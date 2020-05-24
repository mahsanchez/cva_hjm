#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iterator>
#include <iomanip>
#include "product.h"

// (Musiela Parameterisation HJM) We simulate f(t+dt)=f(t) + dfbar   where SDE dfbar =  m(t)*dt+SUM(Vol_i*phi*SQRT(dt))+dF/dtau*dt
void stochastic_process_evolve(const int dimension, const int tenors, double *gaussian_rand, std::vector<std::vector<double>> &fwd_rates, std::vector<double> &drifts, std::vector<std::vector<double>> &volatilities, const int sim, const double dt, double dtau)
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


// Run simulation to generate Forward Rates Risk Factors Grid using HJM Model
void simulate(int dimension, int simN, int tenor_size, double *phi_random, double dt, std::vector<std::vector<double>>& fwd_rates, std::vector<double>& drifts, std::vector<std::vector<double>> volatilities, std::vector<std::vector<double>> exposures, double notional, double K, double maturity, double dtau = 0.5) {
    for (int sim = 1; sim < simN; sim++) {
        // MC path Generation Path Generator or random forward rates risk factors using HJM and Musiela SDE ( CURAND )
        // only simulation points that coincidence with IRS pricing points need to be stored in fwd_rates
        for (int sim = 1; sim < simN; sim++) {
            stochastic_process_evolve(dimension, tenor_size, & phi_random[sim*dimension], fwd_rates, drifts, volatilities, sim, dt, dtau);
        }

        // Vainilla Interest Rate Swap Mark to Market Exposure profile
        InterestRateSwapExposureEngine(fwd_rates, timepoints, notional, K, maturity, dtau, dt).calculate(exposures[sim]);
        //display_curve(exposures[sim]);
    }
}

/*
 *  Reduction to produce the expected exposure profile EE[t] curve for the IRS
 * Expected Exposure  EPE(t) = 𝔼 [max(V , 0)|Ft]
*/
void reduce(std::vector<double>& expected_exposure, std::vector<std::vector<double>>& exposures, int simN, int timepoints_size) {
    for (int t = 0; t < timepoints_size; t++) {
        double sum = 0.0;
        for (int i = 0; i < simN; i += timepoints_size) {
            sum += std::max(exposures[i][t], 0.0);
        }
        expected_exposure[t] = sum /(double)simN;
    }
}


// Monte Carlo Simulation Engine
template<typename GaussianGenerator>
void mc_engine(std::vector<double> &expected_exposure, GaussianGenerator gaussian_generator, int simN, int timepoints_size, double tenor_size,  std::vector<double> &spot_rates, std::vector<std::vector<double>>& volatilities, std::vector<double> &drifts, double expiry = 25.0, double dtau = 0.5) {
    int maturity = 51;
    int dimension = 3;

    double dt = expiry/(double)simN;

    double notional = 1000000; double K = 0.04;

    std::vector<std::vector<double>> exposures(simN, std::vector<double>(tenor_size, 0.0));
    std::vector<std::vector<double>> fwd_rates(simN, std::vector<double>(tenor_size, 0.0));

    // copy the spot rates as the first curve
    std::copy(spot_rates.begin(), spot_rates.end(), fwd_rates[0].begin());

    // Gaussian Random generation
    long randoms_count = dimension * simN;
    double *phi_random = (double *) malloc(sizeof(double) * randoms_count);
    for (int i = 0; i < randoms_count; i++) {
        phi_random[i] = gaussian_generator();
    }

    // Run simulation to generate Forward Rates Risk Factors Grid using HJM Model
    //simulate(dimension, simN, tenor_size, phi_random, dt, fwd_rates, drifts, volatilities, exposures, notional, K, maturity, dtau);

    for (int sim = 1; sim < simN; sim++) {
        // MC path Generation Path Generator or random forward rates risk factors using HJM and Musiela SDE ( CURAND )
        // only simulation points that coincidence with IRS pricing points need to be stored in fwd_rates
        for (int sim = 1; sim < simN; sim++) {
            stochastic_process_evolve(dimension, tenor_size, & phi_random[sim*dimension], fwd_rates, drifts, volatilities, sim, dt, dtau);
        }

        // Vainilla Interest Rate Swap Mark to Market Exposure profile
        InterestRateSwapExposureEngine(fwd_rates, timepoints, notional, K, maturity, dtau, dt).calculate(exposures[sim]);
        //display_curve(exposures[sim]);
    }

    free(phi_random);

    /*
     * Expected Exposure  EPE(t) = 𝔼 [max(V , 0)|Ft]
     * Use reduction to produce the expected exposure profile for the IRS
     */
    // reduce(std::vector<double>& expected_exposure, std::vector<std::vector<double>>& exposures, int simN, int timepoints_size)
    for (int t = 0; t < timepoints_size; t++) {
        double sum = 0.0;
        for (int i = 0; i < simN; i += timepoints_size) {
            sum += std::max(exposures[i][t], 0.0);
        }
        expected_exposure[t] = sum /(double)simN;
    }

    display_curve(expected_exposure);
}


