#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iterator>
#include <iomanip>
#include "product.h"


/*
 *  Reduction to produce the expected exposure profile EE[t] curve for the IRS
 * Expected Exposure  EPE(t) = 𝔼 [max(V , 0)|Ft]
*/
void reduce(std::vector<double>& expected_exposure, std::vector<std::vector<double>>& exposures, int timepoints_size, int simN) {
    for (int t = 0; t < timepoints_size; t++) {
        double sum = 0.0;
        for (int i = 0; i < simN; i++) {
            sum += exposures[i][t];
        }
        expected_exposure[t] = std::max(sum, 0.0) /(double)simN;
    }
}


/*
* Monte Carlo Simulation Engine
* TODO - capture statistics (MC Error, stddeviation, avg)
* Run simulation to generate the stochastics Forward Rates Risk Factors Grid using HJM Model
 */
template<typename GaussianGenerator, typename StochasticProcess>
void mc_engine(std::vector<double> &expected_exposure, GaussianGenerator gaussian_generator, StochasticProcess stochasticProcess, int simN,  double expiry, double dtau) {

    // Simulation Step
    double dt = expiry/(double)simN;

    // Pricing Points Size
    int timepoints_size = expiry/dtau;

    // Interest Rate Swap Pricing Parameters
    double notional = 1000000;
    double K = 0.04;

    // Generate one Exposure Profile by Simulation Step
    std::vector<std::vector<double>> exposures(simN, std::vector<double>(timepoints_size, 0.0));

    // Gaussian Random
    long randoms_count = stochasticProcess.dimension * simN;
    std::vector<double> phi_random(randoms_count, 0.0);

    // Run simulation to generate Forward Rates Risk Factors Grid using HJM Model and Musiela SDE ( CURAND )
    for (int sim = 1; sim < simN; sim++)
    {
        //  Gaussian Random Number Generation
        for (int i = 0; i < randoms_count; i++) {
            phi_random[i] = gaussian_generator();
        }

        // Evolve the Forward Rates Risk Factor Simulation using HJM MonteCarlo
        for (int s = 1; s < simN; s++) {
            for (int t = 0; t < stochasticProcess.tenors_size; t++) {
                stochasticProcess.evolve(s, t, &phi_random[s * stochasticProcess.dimension]);
            }
            // Display Simulation Grid
            // displaySimulationGrid(stochasticProcess.fwd_rates[s], &phi_random[s * stochasticProcess.dimension]);
        }

        // Mark to Market Exposure profile calculation for Interest Rate Swap
        for (int t = 1; t < pricing_points.size() - 1; t++) {
            double reference_day = pricing_points[t];
            HJMYieldCurveTermStructure yieldCurve(stochasticProcess.fwd_rates, reference_day, expiry, dt, dtau);
            VanillaInterestRateSwap irs(notional, t, reference_day, K, pricing_points, yieldCurve);
            exposures[sim][t] = irs.price();
            // Display Discount Factors
            //std::cout << reference_day << "\n"; display_curve(yieldCurve.getDiscountCurve().getPoints());
        }

        // Display EE[t] curve realization for simulation sim
        //display_curve(exposures[sim]);
    }

    /*
     * Expected Exposure  EPE(t) = 𝔼 [max(V , 0)|Ft]
     * Use reduction to produce the expected exposure profile for the IRS
     */
    reduce(expected_exposure, exposures, timepoints_size, simN);

    // Report Expected Exposure Profile Curve
    display_curve(expected_exposure);
}



