#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <functional>
#include <iterator>
#include <iomanip>

#include "product.h"

//#include <boost/math/quadrature/trapezoidal.hpp>


// Random number generation
std::mt19937 mt_rand1(0);
std::mt19937 mt_rand2(100);
std::mt19937 mt_rand3(6000);
std::uniform_real_distribution<> phi1(0.0, 1.0);
std::uniform_real_distribution<> phi2(0.0, 1.0);
std::uniform_real_distribution<> phi3(0.0, 1.0);


void display_curve(std::vector<double> &curve) {
    std::cout << std::setprecision(6)<< std::fixed;
    std::copy(curve.begin(), curve.end(), std::ostream_iterator<double>(std::cout, ","));
    std::cout << std::endl;
}

// We simulate f(t+dt)=f(t) + dfbar   where SDE dfbar =  m(t)*dt+SUM(Vol_i*phi*SQRT(dt))+dF/dtau*dt (Musiela Parameterisation HJM)
auto musiela_sde = [](const double dt, int dimension, const double r0, double drift, std::vector<double> volatilities, std::vector<double> phi, double dF, double dtau)
{
    double mu = drift * dt;
    double vol = 0.0;

    #pragma UNROLL(3)
    for (int i = 0; i < dimension; i++) {
        vol +=  volatilities[i] * phi[i];
    }
    double difussion = vol * std::sqrt(dt);
    double adjustment = (1/dtau) * dF * dt;
    double result = r0 + mu + difussion + adjustment;
    return result;
};

template<typename T>
void stochastic_process(T sde, const int dimension, const int tenors, std::vector<double> &r0, std::vector<double> &r, std::vector<double> &drifts, std::vector<std::vector<double>> &volatilities, std::vector<double> &phi_random, const double dt, double dtau)
{
    double drift;
    double dF;
    std::vector<double> volatility;

    for (int t = 0; t < tenors - 1; t++) {
        drift = drifts[t];

        #pragma UNROLL(3)
        for (int d = 0; d < dimension; d++) {
            volatility.push_back( volatilities[d][t] );
        }

        dF = r0[t+1] - r0[t];
        r[t] = sde(dt, dimension, r0[t], drift, volatility, phi_random, dF, dtau);
    }

    // simulate last tenor
    drift = drifts[tenors-1];
    volatility.clear();
    for (int d = 0; d < dimension; d++) {
        volatility.push_back( volatilities[d][tenors-1] );
    }
    dF = r0[tenors-1] - r0[tenors-2];
    r[tenors-1] = sde(dt, dimension, r0[tenors-1], drift, volatility, phi_random, dF, dtau);
};


// path_generator
// TODO - implement a vector of RandomGenerator variables and MC Variance Reduction
void path_generator(const int pathN, std::vector<std::vector<double>> &fwd_rates, std::vector<double> &spot_rates, std::vector<double> &drifts, std::vector<std::vector<double>> &volatility,  const int tenors,  double dt,  double dtau)
{
    int dimension = 3;
    std::vector<double> rates(tenors, 0.0);
    std::vector<double> fwd(spot_rates);
    std::vector<double> phi(dimension, 0.0);

    // path generator
    for (int sim = 1; sim < pathN; sim++) {
        phi[0] = phi1(mt_rand1);
        phi[1] = phi2(mt_rand2);
        phi[2] = phi3(mt_rand3);

        // evolve the full forward rate curve
        stochastic_process(musiela_sde, dimension, tenors, fwd, rates, drifts, volatility, phi, dt, dtau);

        // store simulated rates in the path
        std::copy(rates.begin(), rates.end(), fwd_rates[sim].begin());

        // reuse the last tenor for next iteration
        std::copy(rates.begin(), rates.end(), fwd.begin());
    }
}


// Run simulations per block, if end of block reached but current_simulations < total_simulations store avg (exposure) on mb_grid[ simulation_block]
// take the expected exposure from the memory block

void mc_engine(std::vector<std::vector<double>>& mm_grid, int simN, int timepoints_size, double tenor_size,  std::vector<double> &spot_rates, std::vector<std::vector<double>>& volatilities, std::vector<double> &drifts, int blockSize = 1000, double dt = 0.01, double dtau = 0.5, bool libor_rate = false) {
    // IRS product information
    double notional = 1000000; // notional
    double K = 0.04; // fixed rates IRS
    int maturity = 51;
    const double tolerance = 1e-6;
    const double LGD = 0.40;

    const int pathN = tenor_size / dt;

    std::vector<std::vector<double>> fwd_rates(pathN, std::vector<double>(timepoints_size, 0.0));
    std::copy(spot_rates.begin(), spot_rates.end(), fwd_rates[0].begin());

    // do Print Monte Carlo Generated Grid
    // gpu - divide the simulation in blocks in such a way that you can increase the number of throughput per threads vs normal strategy and benchmark
    // aplit the simulations on equally blockSizes and use a threadblock to run then TBB store at offset (block_size * n + index) where index = [0..blocksize - 1]
    for (int sim = 0; sim < simN; sim++) {
        // MC path Generation Path Generator (Stochastic Process)
        path_generator(pathN, fwd_rates, spot_rates, drifts, volatilities, timepoints_size, dt, dtau);

        /*
         * IRS Mark to Market and store price in the second grid for [0, maturity], [1, maturity], [2, maturity], ... , [maturity, maturity]
         * move forward in time in fwd[1] and recalculate discount factors TT(exp(-rt))
         * move over timepoints to build the simulation grid
         * Uses vector iterator fwd_rate[i].begin(), fwd_rate[i].end() inside the term structures so you can move among the years
         */
        // TODO - review the Interest Rate Swap pricing and introduce dates

        // IRS pricing with HJM simulation
        int delta = 0;
        for (int cashflow = 0; cashflow < timepoints_size - 1; cashflow++)
        {
           delta += 25; // optimize this waste of cpu cycles
           YieldCurveTermStructure forward_rate(fwd_rates[delta]); // use an iterator here
           DiscountFactorsTermStructure zcb(fwd_rates[delta]);
           InterestRateSwap irs(notional, cashflow + 1, maturity, K, floating_schedule, fixed_schedule, forward_rate, zcb); // start
           mm_grid[sim][cashflow] = irs.price();
        }

        // display exposure profile for IRS on simulation [sim]
        //display_curve(mm_grid[sim]);
    }

}


