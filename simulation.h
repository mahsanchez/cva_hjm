#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iterator>
#include <iomanip>
#include "product.h"

/*
* Monte Carlo Simulation Engine
* MC Apply Variance Reduction
* MC capture simulation statistics (MC Error, stddeviation, avg)
* Run simulation to generate the stochastics Forward Rates Risk Factors Grid using HJM Model
 */
template<typename GaussianGenerator, typename StochasticProcess>
void mc_engine(std::vector<std::vector<double>>  &exposures, GaussianGenerator gaussian_generator, StochasticProcess stochasticProcess, int simN,  double expiry, double dtau, double& duration) {

    // Simulation Step
    double dt = expiry/(double)simN;

    // Pricing Points Size
    int timepoints_size = expiry/dtau + 1;

    // Interest Rate Swap Pricing Parameters
    double notional = 1000000;
    double K = 0.04;

    auto start = std::chrono::high_resolution_clock::now();

    // Gaussian Random
    long randoms_count = stochasticProcess.dimension * simN;
    std::vector<double> phi_random(randoms_count, 0.0);

    // Run simulation to generate Forward Rates Risk Factors Grid using HJM Model and Musiela SDE ( CURAND )
    for (int sim = 1; sim < simN; sim++)
    {
        //  Gaussian Random Number Generation
        gaussian_generator(&phi_random[0], randoms_count);

        // Evolve the Forward Rates Risk Factor Simulation using HJM MonteCarlo
        for (int s = 1; s < simN; s++) {
            for (int t = 0; t < stochasticProcess.tenors_size; t++) {
                stochasticProcess.evolve(s, t, &phi_random[s * stochasticProcess.dimension]);
            }
        }
        //displaySimulationGrid(stochasticProcess.fwd_rates, phi_random);

        // Mark to Market Exposure profile calculation for Interest Rate Swap
        for (int t = 1; t < pricing_points.size() - 1; t++) {
            double reference_day = pricing_points[t];
            HJMYieldCurveTermStructure yieldCurve(stochasticProcess.fwd_rates, reference_day, expiry, dt, dtau);
            VanillaInterestRateSwap irs(notional, t, reference_day, K, pricing_points, yieldCurve);
            exposures[sim][t] = irs.price();
            //displayDiscountFactors(yieldCurve.getDiscountCurve().getPoints());
        }

        // Display EE[t] curve realization for simulation sim
        //display_curve(exposures[sim]);
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    duration = elapsed.count();
}



