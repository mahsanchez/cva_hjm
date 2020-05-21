#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <functional>
#include <iterator>
#include <iomanip>

#include "product.h"

// The random factor in the SDE requires gaussian random values . Hence error function approximate to inverse CDF erfInv(uniform)
// Note : Use allways the same seed value
std::random_device rd{};
std::mt19937 uniform_distribution{rd()};
std::normal_distribution<> gaussian_distribution{0.0,1.0};

// We simulate f(t+dt)=f(t) + dfbar   where SDE dfbar =  m(t)*dt+SUM(Vol_i*phi*SQRT(dt))+dF/dtau*dt (Musiela Parameterisation HJM)
void stochastic_process(const int dimension, const int tenors, std::vector<std::vector<double>> &fwd_rates, std::vector<double> &drifts, std::vector<std::vector<double>> &volatilities, std::vector<double> &phi_random, const int sim, const double dt, double dtau)
{
    double dF;
    double sum_vol;

    for (int t = 0; t < tenors - 1; t++) {
        sum_vol = 0.0;
        for (int i = 0; i < dimension; i++) {
            sum_vol +=  volatilities[i][t] * phi_random[i];
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
        sum_vol +=  volatilities[i][t] * phi_random[i];
    }
    dF = fwd_rates[sim-1][t] - fwd_rates[sim-1][t-1];
    double sde = fwd_rates[sim-1][t];
    sde += drifts[t] * dt;
    sde += sum_vol * std::sqrt(dt);
    sde += (dF/dtau) * dt;
    fwd_rates[sim][t] = sde;
};


// path_generator
void path_generator(const int pathN, std::vector<std::vector<double>> &fwd_rates, std::vector<double> &spot_rates, std::vector<double> &drifts, std::vector<std::vector<double>> &volatility,  const int tenors,  double dt,  double dtau)
{
    int dimension = 3;

    // generate randoms across all paths for SDE
    std::vector<std::vector<double>> phi_randoms(pathN, std::vector<double>(dimension, 0.0));
    for (int i = 0; i < phi_randoms.size(); i++) {
        for(int j = 0; j < dimension; j++) {
            phi_randoms[i][j] = gaussian_distribution(uniform_distribution);
        }
    }

    // path generator evolve the full forward rate curve
    for (int sim = 1; sim < pathN; sim++) {
        stochastic_process(dimension, tenors, fwd_rates, drifts, volatility, phi_randoms[sim], sim, dt, dtau);
    }
}

// Monte Carlo Simulation Engine
void mc_engine(std::vector<double> &expected_exposure, int simN, int timepoints_size, double tenor_size,  std::vector<double> &spot_rates, std::vector<std::vector<double>>& volatilities, std::vector<double> &drifts, int blockSize = 1000, double dt = 0.01, double dtau = 0.5) {
    int maturity = 51;

    const int pathN = tenor_size / dt;

    double notional = 1000000; double K = 0.04;

    std::vector<std::vector<double>> exposures(simN, std::vector<double>(timepoints_size, 0.0));   // TODO - dynamic memory allocation
    std::vector<std::vector<double>> fwd_rates(pathN, std::vector<double>(timepoints_size, 0.0));  // TODO - dynamic memory allocation

    std::copy(spot_rates.begin(), spot_rates.end(), fwd_rates[0].begin());

    // mc_simulation_cpu  - kernel & device code distribute equally block simulations per thread block
    // TODO - Loop Tiling
    for (int sim = 0; sim < simN; sim++) {

        // MC path Generation Path Generator or random forward rates using HJM and Musiela SDE ( CURAND )
        path_generator(pathN, fwd_rates, spot_rates, drifts, volatilities, timepoints_size, dt, dtau);

        // TODO - review the Interest Rate Swap pricing and introduce dates
        // Vainilla Interest Rate Swap Mark to Marking Exposure profile
        InterestRateSwapExposureEngine(exposures[sim], fwd_rates, timepoints, notional, K, maturity, dtau, dt).calculate();
    }

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


