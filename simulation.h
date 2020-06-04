#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iterator>
#include <iomanip>
#include <random>
#include <boost/math/special_functions/erf.hpp>

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

    HJMStochasticProcess(std::vector<double> &spot_rates_, std::vector<double> &drifts_, std::vector<std::vector<double>>& volatilities_, std::vector<std::vector<double>>& fwd_rates_, int dimension_, double dt_, double dtau_, int tenors_size_, double expiry_) :
        spot_rates(spot_rates_), drifts(drifts_), volatilities(volatilities_), fwd_rates(fwd_rates_), dimension(dimension_), dt(dt_), dtau(dtau_), tenors_size(tenors_size_), expiry(expiry_) {

    }

    inline void evolve(int sim, int tenor, double *gaussian_rand) {
        double dfbar = 0.0;

        // calculate diffusion term SUM(Vol_i*phi*SQRT(dt))
        dfbar +=  volatilities[0][tenor] * gaussian_rand[0];
        dfbar +=  volatilities[1][tenor] * gaussian_rand[1];
        dfbar +=  volatilities[2][tenor] * gaussian_rand[2];
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
    std::vector<double> &spot_rates;
    std::vector<std::vector<double>>& volatilities;
    std::vector<std::vector<double>>& fwd_rates;
    int dimension;
    double dt;
    double dtau;
    int tenors_size;
    double expiry;
};

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
    MonteCarloSimulation(PayOff &payOff_, GaussianGenerator &gaussian_generator_, StochasticProcess &stochasticProcess_, std::vector<double>& phi_random_, int simN_) :
            payOff(payOff_), phi_random(phi_random_), gaussian_generator(gaussian_generator_), stochasticProcess(stochasticProcess_), simN(simN_)
     {
        // Simulation Step
        dt = payOff.expiry/(double)simN;
        // Gaussian Random Variates
        randoms_count = simN * stochasticProcess.dimension;
    }

    /**
     * Monte Carlo Method calculation method engine
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
#endif

       // Display EE[t] curve realization for simulation sim
       //display_curve(exposures[sim]);
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        duration = elapsed.count();
    }

    /*
     * Trivial implementation to Evolve the Forward Rates Risk Factor Simulation using HJM MonteCarlo
     */
    void generatePaths() {
        for (int s = 1; s < simN; s++) {
            for (int t = 0; t < stochasticProcess.tenors_size; t++) {
                stochasticProcess.evolve(s, t, &phi_random[s * stochasticProcess.dimension]);
            }
#if DEBUG
            display_curve(stochasticProcess.fwd_rates[s]);
#endif
        }
    }

    /*
     * Mark to Market Exposure profile calculation for Interest Rate Swap (marktomarket) pricePayOff()
     */
    void pricePayOff(std::vector<double> &exposure) {
        HJMYieldCurveTermStructure yieldCurve(stochasticProcess.fwd_rates, payOff.expiry, dt, payOff.dtau);
        double reference_day = payOff.pricing_points[1];
        int pricing_point = 0;
        for (double reference_day = payOff.pricing_points[1]; reference_day < payOff.expiry; reference_day += payOff.dtau) {
            VanillaInterestRateSwapPricer pricer(reference_day, payOff, yieldCurve);
            exposure[pricing_point] = pricer.price();
            pricing_point++;
        }
#if DEBUG
        display_curve(exposure, "IRS Exposure");
#endif
    }

protected:
#if defined(__CUDA_ARCH__)
    curandState * __restrict rngStates;
    double *fwd_rates;
#endif
    std::vector<double>& phi_random;
    PayOff &payOff;
    GaussianGenerator &gaussian_generator;
    StochasticProcess &stochasticProcess;
    int randoms_count;
    int simN;
    double dt;
    std::vector<double> discounts;
    std::vector<double> fra;
};



