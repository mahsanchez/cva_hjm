#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <functional>
#include <iterator>
#include <iomanip>

#include <boost/math/special_functions/erf.hpp>
#include "product.h"

// The random factor in the SDE requires gaussian random values . Hence error function approximate to inverse CDF erfInv(uniform)
// Note : Use allways the same seed value
std::random_device rd{};
/*
std::mt19937 mt{rd()};
std::uniform_real_distribution<double> uniform_distribution{0.0,1.0};
std::normal_distribution<double> gaussian_distribution{0.0,1.0};
*/

// We simulate f(t+dt)=f(t) + dfbar   where SDE dfbar =  m(t)*dt+SUM(Vol_i*phi*SQRT(dt))+dF/dtau*dt (Musiela Parameterisation HJM)
void stochastic_process(const int dimension, const int tenors, double *gaussian_rand, std::vector<std::vector<double>> &fwd_rates, std::vector<double> &drifts, std::vector<std::vector<double>> &volatilities, const int sim, const double dt, double dtau)
{
    double dF;
    double sum_vol;

    for (int t = 0; t < tenors - 1; t++) {
        sum_vol = 0.0;
        for (int i = 0; i < dimension; i++) {
            sum_vol +=  volatilities[i][t] * gaussian_rand[i];
        }
        dF = fwd_rates[sim-1][t+1] - fwd_rates[sim-1][t];
        double sde = fwd_rates[sim-1][t];
        sde += drifts[t] * dt;
        sde += sum_vol * std::sqrt(dt);
        sde += (dF/dtau) * dt;
        fwd_rates[sim][t] = sde;
    }

    // simulate last tenor
    int t = tenors - 1;
    sum_vol = 0.0;
    for (int i = 0; i < dimension; i++) {
        sum_vol +=  volatilities[i][t] * gaussian_rand[i];
    }
    dF = fwd_rates[sim-1][t] - fwd_rates[sim-1][t-1];
    double sde = fwd_rates[sim-1][t];
    sde += drifts[t] * dt;
    sde += sum_vol * std::sqrt(dt);
    sde += (dF/dtau) * dt;
    fwd_rates[sim][t] = sde;
};


// Monte Carlo Simulation Engine
void mc_engine(std::vector<double> &expected_exposure, int simN, int timepoints_size, double tenor_size,  std::vector<double> &spot_rates, std::vector<std::vector<double>>& volatilities, std::vector<double> &drifts, int blockSize = 1000, double dt = 0.01, double dtau = 0.5) {
    int maturity = 51;
    int dimension = 3;

    const int pathN = tenor_size / dt;

    double notional = 1000000; double K = 0.04;

    std::vector<std::vector<double>> exposures(simN, std::vector<double>(timepoints_size, 0.0));
    std::vector<std::vector<double>> fwd_rates(pathN, std::vector<double>(timepoints_size, 0.0));

    // copy the spot rates as the first curve
    std::copy(spot_rates.begin(), spot_rates.end(), fwd_rates[0].begin());

    auto real_rand = std::bind(std::uniform_real_distribution<double>(0.0,1.0), mt19937(rd()));

    // Gaussian Random generation
    double *phi_random = (double *) malloc(sizeof(double) * dimension * pathN);
    for (int i = 0; i < dimension * pathN; i++) {
        phi_random[i] = boost::math::erf_inv(real_rand());
    }

    // mc_simulation_cpu  - kernel & device code distribute equally block simulations per thread block
    // TODO - Loop Tiling
    for (int sim = 1; sim < simN; sim++) {

        // MC path Generation Path Generator or random forward rates risk factors using HJM and Musiela SDE ( CURAND )
        for (int sim = 1; sim < pathN; sim++) {
            stochastic_process(dimension, tenor_size, & phi_random[sim*dimension], fwd_rates, drifts, volatilities, sim, dt, dtau);
        }

        // Vainilla Interest Rate Swap Mark to Marking Exposure profile
        InterestRateSwapExposureEngine(exposures[sim], fwd_rates, timepoints, notional, K, maturity, dtau, dt).calculate();
        //display_curve(exposures[sim]);
    }

    free(phi_random);

    /*
     * Expected Exposure  EPE(t) = 𝔼 [max(V , 0)|Ft]
     * TODO - USe parallel reduction to produce the expected exposure profile for the IRS
     */
    for (int t = 0; t < timepoints_size; t++) {
        double sum = 0.0;
        for (int i = 0; i < simN; i += timepoints_size) {
            sum += std::max(exposures[i][t], 0.0);
        }
        expected_exposure[t] = sum /(double)simN;
    }

    display_curve(expected_exposure);
}


