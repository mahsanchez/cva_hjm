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


// Run simulation to generate the stochastics Forward Rates Risk Factors Grid using HJM Model
// TODO - capture statistics (MC Error, stddeviation, avg)
void simulate(int dimension, int simN, int tenor_size, double *phi_random, double dt, std::vector<std::vector<double>>& fwd_rates, std::vector<double>& drifts, std::vector<std::vector<double>> volatilities, std::vector<std::vector<double>>& exposures, double notional, double K, double expiry, double dtau = 0.5) {
    for (int sim = 1; sim < simN; sim++) {
        // MC path Generation Path Generator or random forward rates risk factors using HJM and Musiela SDE ( CURAND )
        // only simulation points that coincidence with IRS pricing points need to be stored in fwd_rates
        for (int s = 1; s < simN; s++) {
            stochastic_process_evolve(dimension, tenor_size, & phi_random[sim*dimension], fwd_rates, drifts, volatilities, s, dt, dtau);
        }

        // Vainilla Interest Rate Swap Mark to Market Exposure profile (Interest Rate Swap)
        InterestRateSwapExposureEngine(fwd_rates, timepoints, notional, K, expiry, dtau, dt).calculate(exposures[sim]);

        //display_curve(exposures[sim]);
    }
}

/*
 *  Reduction to produce the expected exposure profile EE[t] curve for the IRS
 * Expected Exposure  EPE(t) = 𝔼 [max(V , 0)|Ft]
*/
void reduce(std::vector<double>& expected_exposure, std::vector<std::vector<double>>& exposures, int timepoints_size, int simN) {
    for (int t = 0; t < timepoints_size; t++) {
        double sum = 0.0;
        for (int i = 0; i < simN; i += timepoints_size) {
            sum += std::max(exposures[i][t], 0.0);
        }
        expected_exposure[t] = sum /(double)simN;
    }
}


// Monte Carlo Simulation Engine
template<typename GaussianGenerator> //<typename GaussianGenerator, StochasticProcess>
void mc_engine(std::vector<double> &expected_exposure, GaussianGenerator gaussian_generator, int simN, int timepoints_size, std::vector<double> &spot_rates, std::vector<std::vector<double>>& volatilities, std::vector<double> &drifts, double expiry = 25.0, double dtau = 0.5) {
    int dimension = 3;

    double dt = expiry/(double)simN;

    double notional = 1000000; double K = 0.04;

    // Generate a grid of tenor_size x simN contains the risk factors Forward Rates F(t, T1, T2) from 0 .. T (expiry)
    std::vector<std::vector<double>> fwd_rates(simN, std::vector<double>(timepoints_size, 0.0));

    // Generate one Exposure Profile by Simulation Step
    std::vector<std::vector<double>> exposures(simN, std::vector<double>(timepoints_size, 0.0));

    // copy the spot rates as the first curve
    std::copy(spot_rates.begin(), spot_rates.end(), fwd_rates[0].begin());

    // Gaussian Random generation
    long randoms_count = dimension * simN;
    double *phi_random = (double *) malloc(sizeof(double) * randoms_count);
    for (int i = 0; i < randoms_count; i++) {
        phi_random[i] = gaussian_generator();
    }

    // Run simulation to generate Forward Rates Risk Factors Grid using HJM Model
    simulate(dimension, simN, timepoints_size, phi_random, dt, fwd_rates, drifts, volatilities, exposures, notional, K, expiry, dtau);

    // Release resources
    free(phi_random);

    /*
     * Expected Exposure  EPE(t) = 𝔼 [max(V , 0)|Ft]
     * Use reduction to produce the expected exposure profile for the IRS
     */
    reduce(expected_exposure, exposures, timepoints_size, simN);

    // Report Expected Exposure Profile Curve
    display_curve(expected_exposure);
}


