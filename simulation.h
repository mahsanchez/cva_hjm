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


std::vector<double> timepoints = {0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75,
                            9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25, 11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0, 14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 17.75, 16.0, 16.25, 16.5, 16.75,
                            17.0, 17.25, 17.5, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5, 19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25, 22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25,  24.5, 24.75, 25.0 };


std::vector<double> floating_schedule = {0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75,
                                   9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25, 11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0, 14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 17.75, 16.0, 16.25, 16.5, 16.75,
                                   17.0, 17.25, 17.5, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5, 19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25, 22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0 };

std::vector<double> fixed_schedule = {0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75,
                                         9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25, 11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0, 14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 17.75, 16.0, 16.25, 16.5, 16.75,
                                         17.0, 17.25, 17.5, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5, 19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25, 22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0 };



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
    std::vector<double> phi;

    // path generator
    for (int sim = 1; sim < pathN; sim++) {
        phi.push_back(phi1(mt_rand1));
        phi.push_back(phi2(mt_rand2));
        phi.push_back(phi3(mt_rand3));

        // evolve the full forward rate curve
        stochastic_process(musiela_sde, dimension, tenors, fwd, rates, drifts, volatility, phi, dt, dtau);
        // store simulated rates in the path
        fwd_rates.push_back(rates);
        // reuse the last tenor for next iteration
        std::copy(rates.begin(), rates.end(), fwd.begin());

        //clear randoms
        phi.clear();
    }
}


double mc_engine(int simN, int timepoints_size, double tenor_size,  std::vector<double> &spot_rates, std::vector<std::vector<double>>& volatilities, std::vector<double> &drifts, int blockSize = 1000, double dt = 0.01, double dtau = 0.5, bool libor_rate = false) {
    // IRS product information
    double notional = 1000000; // notional
    double K = 0.04; // fixed rates IRS
    int maturity = 51;
    const double tolerance = 1e-6;
    const double LGD = 0.40;

    double cva = tolerance;

    std::vector<std::vector<double>> fwd_rates;
    std::vector<std::vector<double>> mm_grid;

    // do Print Monte Carlo Generated Grid
    // gpu - divide the simulation in blocks in such a way that you can increase the number of throughput per threads vs normal strategy and benchmark
    // aplit the simulations on equally blockSizes and use a threadblock to run then TBB store at offset (block_size * n + index) where index = [0..blocksize - 1]
    for (int sim = 0; sim < simN; sim++) {
        // MC path Generation Path Generator (Stochastic Process)
        const int pathN = tenor_size / dt;
        fwd_rates.push_back(spot_rates);

        path_generator(pathN, fwd_rates, spot_rates, drifts, volatilities, timepoints_size, dt, dtau);

        /*
         * IRS Mark to Market and store price in the second grid for [0, maturity], [1, maturity], [2, maturity], ... , [maturity, maturity]
         * move forward in time in fwd[1] and recalculate discount factors TT(exp(-rt))
         * move over timepoints to build the simulation grid
         * Uses vector iterator fwd_rate[i].begin(), fwd_rate[i].end() inside the term structures so you can move among the years
         */
        // TODO - review the Interest Rate Swap pricing and introduce dates
        std::vector<double> mark_to_marking(timepoints_size);

        // For each tenors across timepoints
        int delta = 0;
        for (int cashflow = 0; cashflow < timepoints_size - 1; cashflow++)
        {
           delta += 25; // optimize this waste of cpu cycles
           YieldCurveTermStructure forward_rate(fwd_rates[delta]); // use an iterator here
           DiscountFactorsTermStructure zcb(fwd_rates[delta]);
           InterestRateSwap irs(notional, cashflow + 1, maturity, K, floating_schedule, fixed_schedule, forward_rate, zcb); // start
           mark_to_marking[cashflow] = irs.price();
        }

        // add exposure curve to exposure profile grid
        mm_grid.push_back(mark_to_marking);

        // remove fwd rates simulation
        fwd_rates.clear();

        // Display in the stdout the simulated Curve
        //display_curve(mark_to_marking);
    }


    // Run simulations per block, if end of block reached but current_simulations < total_simulations store avg (exposure) on mb_grid[ simulation_block]
    // take the expected exposure from the memory block
    std::vector<double> pd(timepoints_size, 0.01);
    std::vector<double> eexposure(timepoints_size, 0.00);

    // calculate expected exposure profile  Average [ Max ( irs[i], 0) ] , n * 51 reductions
    // gpu - store the Mark to Market curve in shared memory and avg them by block, apply reduction to avg values
    // gpu - store the Mark to Market curve in global memory and avg the value on each datapoints
    const double factor = 1.0/(double)simN;
    for (int t = 0; t < timepoints_size; t++) {
        double sum = 0.0;
        for (int sim = 0; sim < simN; sim++) {
            sum += ( mm_grid[sim][t] > 0) ? mm_grid[sim][t] : 0.0;
        }
        eexposure[t] = sum * factor;
    }

    // Clear the mm_grid simulation for next iteration
    mm_grid.clear();

    // Display in the stdout the simulated forward rate
    //display_curve(eexposure);

    // TODO - implement probability of default calculation (pd)
    // cpu - single thread implementation

    // Generate CVA = LGD * INTEGRAL (EE(t), PD(t), DF(t)) (T0, T) done by lambda Apply Trapezoidal Rule integration
    TermStructure ee_curve( eexposure );
    TermStructure pd_curve( pd );
    DiscountFactorsTermStructure df_curve(spot_rates);

    // integration using trapezoidal curve
    cva = CVAPricer(LGD, eexposure, pd, df_curve).price();

    return cva;
}


