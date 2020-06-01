#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iterator>
#include <iomanip>
#include <random>
//#include <boost/math/special_functions/erf.hpp>

#if defined(__CUDA_ARCH__)
    #include <cuda_runtime_api.h>
#endif

#include "product.h"

#if defined(__CUDA_ARCH__)
#define checkError(code) {
     if ((code) != cudaSuccess) {
        fprintf(stderr, "Cuda failure %s:%d: '%s' \n", __FILE__, __LINE__, cudaGetErrorString(code));
     }
}
#endif

/*
 * Gaussian Random Number generators
 */
std::random_device rd;
std::mt19937  mt(rd() );
std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);
std::normal_distribution<double> gaussian_distribution(0.0, 1.0);

// Mike Giles ErfInv Aproximation of erfInv using a polinomio interpolation
// https://people.maths.ox.ac.uk/gilesm/files/gems_erfinv.pdf available at GPU Gems
#if defined(__CUDA_ARCH__)
__inline__ __device__
#endif
double MBG_erfinv(double x)
{
    float w, p;

#if defined(__CUDA_ARCH__)
    w = - __logf((1.0f-x)*(1.0f+x));
#else
    w = - std::log((1.0f-x)*(1.0f+x));
#endif
    if ( w < 5.000000f ) {
        w = w - 2.500000f;
        p = 2.81022636e-08f;
        p = 3.43273939e-07f + p*w;
        p = -3.5233877e-06f + p*w;
        p = -4.39150654e-06f + p*w;
        p = 0.00021858087f + p*w;
        p = -0.00125372503f + p*w;
        p = -0.00417768164f + p*w;
        p = 0.246640727f + p*w;
        p = 1.50140941f + p*w;
    }
    else {
#if defined(__CUDA_ARCH__)
        w = sqrtf(w) - 3.000000f;
#else
        w = std::sqrt(w) - 3.000000f;
#endif
        p = -0.000200214257f;
        p = 0.000100950558f + p*w;
        p = 0.00134934322f + p*w;
        p = -0.00367342844f + p*w;
        p = 0.00573950773f + p*w;
        p = -0.0076224613f + p*w;
        p = 0.00943887047f + p*w;
        p = 1.00167406f + p*w;
        p = 2.83297682f + p*w;
    }
    return p*x;
}

//Error Function Inverse over uniform distributed numbers to generate Gaussian Random Generators
class ErfInvGaussianRandomGenerator {
public:
    void operator() (double *phi_random, int count) {
        for (int i = 0; i < count; i++) {
            //phi_random[i] = boost::math::erf_inv(uniform_distribution(mt));
            phi_random[i] = MBG_erfinv(uniform_distribution(mt));
        }
    }
};

// Gaussian generator using C++ STL normal_distribution
class NormalRandomGenerator {
public:
    void operator() (double *phi_random, int count) {
        for (int i = 0; i < count; i++) {
            phi_random[i] = gaussian_distribution(mt);
        }
    }
};

/*
 * (Musiela Parameterisation HJM) We simulate f(t+dt)=f(t) + dfbar  where SDE dfbar =  m(t)*dt+SUM(Vol_i*phi*SQRT(dt))+dF/dtau*dt  and phi ~ N(0,1)
 */
struct HJMStochasticProcess {

    HJMStochasticProcess(std::vector<double> &drifts_, std::vector<std::vector<double>>& volatilities_, std::vector<std::vector<double>>& fwd_rates_, int dimension_, double dt_, double dtau_, int tenors_size_, double expiry_) :
            drifts(drifts_), volatilities(volatilities_), fwd_rates(fwd_rates_), dimension(dimension_), dt(dt_), dtau(dtau_), tenors_size(tenors_size_), expiry(expiry_) {

    }

    inline void evolve(int sim, int tenor, double *gaussian_rand) {
        double dfbar = 0.0;

        // calculate diffusion term SUM(Vol_i*phi*SQRT(dt))
        for (int i = 0; i < dimension; i++) {
            dfbar +=  volatilities[i][tenor] * gaussian_rand[i];
        }
        dfbar *= std::sqrt(dt);

        // calculate the drift m(t)*dt
        dfbar += drifts[tenor] * dt;

        // calculate dF/dtau*dt
        double dF = 0.0;
        if (tenor < (tenors_size - 1)) {
            dF += fwd_rates[sim-1][tenor+1] - fwd_rates[sim-1][tenor];
        }
        else {
            dF += fwd_rates[sim-1][tenor] - fwd_rates[sim-1][tenor-1];
        }
        dfbar += (dF/dtau) * dt;

        // apply Euler Maruyana
        dfbar = fwd_rates[sim-1][tenor] + dfbar;

        fwd_rates[sim][tenor] = dfbar;
    }

    std::vector<double> &drifts;
    std::vector<std::vector<double>>& volatilities;
    std::vector<std::vector<double>>& fwd_rates;
    int dimension;
    double dt;
    double dtau;
    int tenors_size;
    double expiry;
};


/*
 * This kernel initialize CURAND gaussian random generator.
 * place each device pathN random numbers apart on the random number sequence
 */
#if defined(__CUDA_ARCH__)
static __global__ void rngSetupStates( curandState *rngState, int device_id)
{
    // determine global thread id
    int tid = threadIdx.x + blockIdx.x * blockDim.x;  //f1pNTu8LT!
    // Each threadblock gets different seed,
    // Threads within a threadblock get different sequence numbers
    curand_init(blockIdx.x + gridDim.x * device_id, threadIdx.x, 0, &rngState[tid]);
}

/*
 * (Musiela Parameterisation HJM) We simulate f(t+dt)=f(t) + dfbar  where SDE dfbar =  m(t)*dt+SUM(Vol_i*phi*SQRT(dt))+dF/dtau*dt  and phi ~ N(0,1)
 * This SDE HJM is evolved like a stencil. Each point value is calculated in a point by resting the value of the upper right point neighbor minus the uppter left neighbor point
 * Latest point value is totally the other way around
 */
__global__
inline double sde_hjm_evolve(int sim, int tenor, double *fwd_rates, double *gaussian_rand, double *volatilities, double *drift) {
        double dfbar = 0.0;
        for (int i = 0; i < dimension; i++) {
            dfbar +=  volatilities[i][tenor] * gaussian_rand[i];
        }
        dfbar *= std::sqrt(dt);
        dfbar += drifts[tenor] * dt;
        double dF = 0.0;
        if (tenor < (tenors_size - 1)) {
            dF += fwd_rates[tenor+1] - fwd_rates[tenor];
        }
        else {
            dF += fwd_rates[tenor] - fwd_rates[tenor-1];
        }
        dfbar += (dF/dtau) * dt;
        dfbar = fwd_rates[sim-1][tenor] + dfbar;
        return dfbar;
}

/*
 * This kernel computes the integral over all paths using a single thread block per simulation block.
 * It is fastest when the number of thread blocks times the work per block is high enough to keep the GPU busy.
 */
__global__
void __generatePaths(double* exposure, double* forward_rates, double *volatilities, double *drifts, curandState * __restrict rngStates, int sim, int pathN, int chunk, double expiry, double dt, double dtau) {
    const int globalIdx = blockIdx.x * blockDim.x + threadIdx.x;
    const int t = threadIdx.x;
    const int points = expiry/dt + 1;
    const int dimension = stochasticProcess.dimension;
    double phi_random[dimension];
    __shared__ double _fwdCurve[points*pathN];

    for (int i = 0; i < chunk_size; i++) {
        // Copy random number state to local memory for efficiency
        curandState localState = rngStates[tid];
        for (int i = 0; i < dimension * pathN; i++) {
            phi_random[i] = curand_normal(&localState);
        }

        //Initialize Initial fwdCurve Shared Memory
        _fwdCurve[ threadIdx.x ] = forward_rates[blockIdx.x * blockDim.x + threadIdx.x];
        __syncthreads();

        //Cycle through the entire samples array forward rate for each path accumulate the rates values into intermediate shared memory buffer
        for (int s = 1; s < pathN; s++) {
            _fwdCurve[s * pathN + threadIdx.x] = sde_hjm_evolve(s, threadIdx.x, phi_random, volatilities, drift);
            __syncthreads();
        }

        // Cycle through the entire random paths array and copy the fwdCurve to Global Memory
        for (int s = 1; s < pathN; s++) {
            forward_rates[s][threadIdx.x] = _fwdCurve[s * pathN + threadIdx.x];
        }

        // Store the partial results over the partial_values with reduction
        // Perform prefixsum to generate the values for the YieldCurve.Discounts
        // Register for interesting simulation point to produce the ForwardRates
    }
}

/*
 * Interest Rate Swap Mark to Market pricing the IRS across all pricing points. Each Thread produce a point for the Exposure Curve
 * This kernel computes the integral over all paths using a single thread block per simulation block.
 * It is fastest when the number of thread blocks times the work per block is high enough to keep the GPU busy.
 */
__global__
void __pricePayOff(double* exposure, double* forward_rates, int sim, int pathN, int chunk_size, double expiry, double dt, double dtau) {
    const int globalIdx = blockIdx.x * blockDim.x + threadIdx.x;
    const int t = threadIdx.x;
    const int points = expiry/dt + 1;
    __shared__ double _fwdCurve[points*pathN];

    for (int i = 0; i < chunk_size; i++) {
        if (threadIdx.x < payOff.expiry/payOff.dtau) {
            double reference_day = payOff.pricing_points[threadIdx.x];
            HJMYieldCurveTermStructure yieldCurve(fwd_rates, reference_day, payOff.expiry, dt, payOff.dtau);
            VanillaInterestRateSwapPricer pricer(t, reference_day, payOff, yieldCurve);
            exposure[i][threadIdx.x] = pricer.price();
        }
    }
}
#endif


/*
* Monte Carlo Simulation Engine
* Run simulation to generate the stochastics Forward Rates Risk Factors Grid using HJM Model
* Simulates forward rates using Euler-Maruyama time-stepping procedure.
* MC Apply Variance Reduction
* MC capture simulation statistics (MC Error, stddeviation, avg)
 */
template<typename GaussianGenerator, typename StochasticProcess, typename PayOff> //Statistics
class MonteCarloSimulation {
public:
    MonteCarloSimulation(PayOff &payOff_, GaussianGenerator &gaussian_generator_, StochasticProcess &stochasticProcess_, int simN_) :
            payOff(payOff_), gaussian_generator(gaussian_generator_), stochasticProcess(stochasticProcess_), simN(simN_)
     {
        // Simulation Step
        dt = payOff.expiry/(double)simN;

        // Gaussian Random values container
        randoms_count = stochasticProcess.dimension * simN;
        phi_random.resize(randoms_count);
    }

    /**
     * Monte Carlo Method calculation method engine
     * Introducing CUDA Streams to alternate between computation and data
     */
    void calculate(std::vector<std::vector<double>> &exposures, double &duration) {
        auto start = std::chrono::high_resolution_clock::now();

        // Run simulation to generate Forward Rates Risk Factors Grid using HJM Model and Musiela SDE ( CURAND )
        #if !defined(__CUDA_ARCH__)
            for (int sim = 1; sim < simN; sim++) {
                // Gaussian Random Number Generation
                gaussian_generator(&phi_random[0], randoms_count);
                // Evolve the Forward Rates Risk Factor Simulation Path using HJM Model
                generatePaths();
                // Interest Rate Swap Mark to Market pricing the IRS across all pricing points
                pricePayOff(exposures[sim]);
            }
        #else
            // Dimension the cuda kernel
            int blockSize = 32;  // Number of Threads per block depends on the number of pricing points for the IRS (irs.expiry/dtau + 1)
            int numBlocks = (simN/pathN + blockSize - 1) / blockSize;
            int SM = 10; // Get this value from the CUDA Device properties, also check the number of active blocks in the device
            int pathN = simN/SM; // Number of concurrent simulation running at a time across all SM

            // Memory input data to be used inside the CUDA device
            double *spot_rates_device;
            double *exposure_device;

            // Allocate device memory to copy the spot rates
            checkCudaErrors(cudaMalloc(&spot_rates_device, tenor_size * sizeof(double));;
            checkCudaErrors(cudaMalloc(&exposure_device, pathN * (payOff.expiry/payOff.dtau + 1) * sizeof(double));;
            checkCudaErrors(cudaMalloc(&forward_rates, pathN * (payOff.expiry/payOff.dtau + 1) * sizeof(double));;
            checkCudaErrors(cudaMalloc((void **) rngStates, numBlocks * blockSize * sizeof(curandState)));

            // Initialise RNG
            rngSetupStates<<<numBlocks, blockSize>>>(rngStates, device);

            // Generate Paths instantaneus forward rate simulation
            __generatePaths<<<numBlocks, blockSize>>>(forward_rates, spot_rates_device, stochasticProcess.volatilities, stochasticProcess.drifts, rngStates, sde_hjm_evolve, sim, pathN, expiry, dt, dtau);

            // Interest Rate Swap Mark To Market and produce the exposure grids
            __pricePayOff<<<numBlocks, blockSize>>>(exposures, forward_rates, sim, pathN, chunk_size, expiry, dt, dtau);

            // Copy the results exposures back to host
            checkCudaErrors(cudaMemcpy(exposures, exposures_device, bytes, cudaMemcpyDeviceToHost));

            // Wait for all gpu jobs to finish
            checkCudaErrors(cudaDeviceSynchronize());

            // CleanUp
            checkCudaErrors(cudaFree(spot_rates));
            checkCudaErrors(cudaFree(forward_rates));
            checkCudaErrors(cudaFree(rngStates));
       #endif

       // Display EE[t] curve realization for simulation sim
       //display_curve(exposures[sim]);

        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        duration = elapsed.count();
    }

    /*
     * Evolve the Forward Rates Risk Factor Simulation using HJM MonteCarlo
     */
    void generatePaths() {
        for (int s = 1; s < simN; s++) {
            for (int t = 0; t < stochasticProcess.tenors_size; t++) {
                stochasticProcess.evolve(s, t, &phi_random[s * stochasticProcess.dimension]);
            }
        }
        //displaySimulationGrid(stochasticProcess.fwd_rates, phi_random);
    }

    /*
     * Mark to Market Exposure profile calculation for Interest Rate Swap (marktomarket) pricePayOff()
     */
    void pricePayOff(std::vector<double> &exposure) {
        for (int t = 1; t < payOff.pricing_points.size() - 1; t++) {
            double reference_day = payOff.pricing_points[t];
            HJMYieldCurveTermStructure yieldCurve(stochasticProcess.fwd_rates, reference_day, payOff.expiry, dt, payOff.dtau);
            VanillaInterestRateSwapPricer pricer(t, reference_day, payOff, yieldCurve);
            exposure[t] = pricer.price();
            //displayDiscountFactors(yieldCurve.getDiscountCurve().getPoints());
        }
    }

protected:
#if defined(__CUDA_ARCH__)
    curandState * __restrict rngStates;
    double *fwd_rates;
#endif
    std::vector<double> phi_random;
    PayOff &payOff;
    GaussianGenerator &gaussian_generator;
    StochasticProcess &stochasticProcess;
    int randoms_count;
    int simN;
    double dt;
};



