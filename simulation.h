#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iterator>
#include <iomanip>
#include <random>
#include <boost/math/special_functions/erf.hpp>

#if defined(__CUDA_ARCH__)
#include <cuda.h>
#endif


#include "product.h"

/*
 * Gaussian Random Number generators
 */
std::random_device rd;
std::mt19937  mt(rd() );
std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);
std::normal_distribution<double> gaussian_distribution(0.0, 1.0);

//Error Function Inverse over uniform distributed numbers to generate Gaussian Random Generators
class ErfInvGaussianRandomGenerator {
public:
    void operator() (double *phi_random, int count) {
        for (int i = 0; i < count; i++) {
            phi_random[i] = boost::math::erf_inv(uniform_distribution(mt));
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
void _simulate_kernel(double *forward_rates, double *volatilities, double *drifts, curandState * __restrict rngStates, int sim, int pathN, double expiry, double dt, double dtau) {
    const int globalIdx = blockIdx.x * blockDim.x + threadIdx.x;
    const int x = threadIdx.x;
    const int points = expiry/dt + 1;
    const int dimension = stochasticProcess.dimension;
    double phi_random[dimension];
    __shared__ double _fwdCurve[points*pathN];

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

        #if defined(__CUDA_ARCH__)
        cudaMalloc(&fwd_rates, pathN * (payOff.expiry/payOff.dtau + 1) * sizeof(double));fwd_rates;
        #endif
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
                //  Gaussian Random Number Generation
                gaussian_generator(&phi_random[0], randoms_count);
                // Evolve the Forward Rates Risk Factor Simulation Path using HJM Model
                simulate();
                // Interest Rate Swap Mark to Market pricing the IRS across all pricing points
                pricePayOff(exposures[sim]);
            }
        #else
            // Dimension the cuda kernel
            int blockSize = 32;  // Number of Threads per block depends on the number of pricing points for the IRS (irs.expiry/dtau + 1)
            int numBlocks = (simN/pathN + blockSize - 1) / blockSize;

            //Allocate states for pseudo random number generators
             checkCudaErrors(cudaMalloc((void **) rngStates, numBlocks * blockSize * sizeof(curandState)));

            // place each device pathN random numbers apart on the random number sequence
            rngSetupStates<<<numBlocks, blockSize>>>(rngStates, device);

            for (int sim = 1; sim < simN; sim++) {
                // launch the cuda kernel
                for (int i = 0; i < simN/pathN; i++) {
                    _simulate_kernel<<<numBlocks, blockSize>>>(forward_rates, stochasticProcess.volatilities, stochasticProcess.drifts, rngStates, sde_hjm_evolve, sim, pathN, expiry, dt, dtau);
                    cudaMemcpy(stochasticProcess.fwd_rates, fwd_rates[index], bytes, cudaMemcpyDeviceToHost);
                }
                // Wait for GPU to finish before accessing on host
                cudaDeviceSynchronize();

                // Price PayOff
                pricePayOff(exposures[sim]);
            }

            // CleanUp
            cudaFree(fwd_rates);
            cudaFree(rngStates));
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
    void simulate() {
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
    double **fwd_rates;
    int pathN = 1024;    // Number of Simulations Blocks must be multiple of 1024
#endif
    std::vector<double> phi_random;
    PayOff &payOff;
    GaussianGenerator &gaussian_generator;
    StochasticProcess &stochasticProcess;
    int randoms_count;
    int simN;
    double dt;
};



