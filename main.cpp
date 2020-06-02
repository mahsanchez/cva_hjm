#include <iostream>
#include <vector>
#include <iterator>

//#include "calibration.h"
#include "simulation.h"

using namespace std;

const int tenor_size = 51;

// First row is the last observed forward curve (BOE data)  //  3-dimensional Normal random vector in columns BC, BD, BE (on far right)
std::vector<double> spot_rates = {0.046138361,0.045251174,0.042915805,0.04283311,0.043497719,0.044053792,0.044439518,0.044708496,0.04490347,0.045056615,0.045184474,0.045294052,0.045386152,0.045458337,0.045507803,0.045534188,0.045541867,0.045534237,0.045513128,0.045477583,0.04542292,0.045344477,0.04523777,0.045097856,0.044925591,0.04472353,0.044494505,0.044242804,0.043973184,0.043690404,0.043399223,0.043104398,0.042810688,0.042522852,0.042244909,0.041978295,0.041723875,0.041482518,0.04125509,0.041042459,0.040845492,0.040665047,0.040501255,0.040353009,0.040219084,0.040098253,0.039989288,0.039890964,0.039802053,0.039721437,0.03964844};

// Volatility Calibration
std::vector<std::vector<double>> volatilities = //(3, std::vector<double>(51, 0.0));
{
{0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655,0.006430655 },
{-0.003556543,-0.003811644,-0.004010345,-0.004155341,-0.004249328,-0.004295001,-0.004295055,-0.004252185,-0.004169088,-0.004048459,-0.003892993,-0.003705385,-0.003488331,-0.003244526,-0.002976667,-0.002687447,-0.002379563,-0.00205571,-0.001718584,-0.001370879,-0.001015292,-0.000654518,-0.000291251,7.18113E-05,0.000431975,0.000786544,0.001132823,0.001468117,0.001789731,0.002094968,0.002381133,0.002645532,0.002885468,0.003098246,0.003281171,0.003431548,0.00354668,0.003623873,0.00366043,0.003653657,0.003600859,0.003499339,0.003346403,0.003139354,0.002875498,0.002552139,0.002166581,0.00171613,0.00119809,0.000609765,-5.15406E-05 },
{-0.004750672,-0.003908573,-0.003134891,-0.00242728,-0.001783395,-0.001200891,-0.000677421,-0.00021064,0.000201797,0.000562236,0.000873023,0.001136502,0.00135502,0.001530923,0.001666555,0.001764262,0.00182639,0.001855285,0.001853291,0.001822755,0.001766022,0.001685437,0.001583346,0.001462095,0.00132403,0.001171495,0.001006836,0.000832399,0.00065053,0.000463573,0.000273875,8.37816E-05,-0.000104363,-0.000288212,-0.00046542,-0.000633643,-0.000790533,-0.000933746,-0.001060936,-0.001169758,-0.001257866,-0.001322914,-0.001362557,-0.001374449,-0.001356245,-0.0013056,-0.001220167,-0.001097601,-0.000935557,-0.000731689,-0.000483652}
};

// Drift Callibration
std::vector<double> drifts = {0.000000,0.000036,0.000069,0.000099,0.000128,0.000155,0.000182,0.000208,0.000233,0.000257,0.000280,0.000303,0.000324,0.000345,0.000364,0.000381,0.000397,0.000412,0.000426,
 0.000438,0.000449,0.000459,0.000469,0.000478,0.000487,0.000496,0.000505,0.000515,0.000526,0.000538,0.000551,0.000566,0.000583,0.000601,0.000621,0.000643,0.000667,0.000692,0.000718,0.000746,0.000774,
 0.000803,0.000832,0.000861,0.000889,0.000916,0.000943,0.000967,0.000991,0.001013,0.001035
};

// Spreads
std::vector<double> spreads = {
        208.5, 187.5, 214, 235, 255, 272, 225, 233, 218, 215, 203, 199, 202, 196.71, 225.92, 219.69, 229.74, 232.12, 223.02, 224.45, 212, 211.51, 206.25, 203.37, 212.94, 211.02, 210.06, 206.23, 211.49, 209.09, 204.3, 204.77,
        199.98, 200.94, 202.38, 205.72, 204.76, 210.02, 209.54, 209.54, 212.41, 213.35, 208.57, 208.56, 220.05, 226.26, 227.2, 222.89, 228.63, 231.5, 247.75, 255.37, 251.07,
        256.33, 252.01, 254.88, 246.98, 238.12, 241.95, 244.33, 252.45, 250.53, 246.71, 256.26, 255.78, 257.19, 247.63, 237.12, 234.73, 226.36, 218, 219.9, 224.68,
        221.81, 220.38, 211.77, 203.17, 206.04, 220, 218, 225, 217.5, 215, 220.5, 250.25, 260.5, 269.5, 258, 258.5, 263, 274.5, 265.5, 268.5, 273.5, 275, 271.5, 263.75, 275, 287, 281, 271, 280.25, 284.5, 272, 275.5, 264.25, 274.25,
        269.5, 264.5, 256.5, 258, 265, 260, 262.5, 268, 272.25, 271, 274, 278.5, 278.5, 284.5, 290, 274, 264.5, 262.5, 247.75, 250.5, 248, 244.75, 246.5, 237.5, 240.5, 236,
        245.5, 237.75, 234.25, 235, 224, 215.5, 217, 217, 220.5, 208.5, 202.47, 203.43, 206.77, 205.33, 200.05, 202.91, 205.05, 222.51, 218.9, 218.43, 221.31, 217.24, 218.67,
        216.52, 216.5, 217.94, 208.37, 205.01, 200.95, 203.1, 203.81, 206.2, 204.28, 200.93, 202.36, 200.44, 197.8, 199.23, 209.74, 211.18, 214.05, 215, 228.62, 233.63, 222.86, 218.8, 214.49, 217.36, 216.4, 213.52, 215.43,
        219.49, 214.22, 218.05, 211.1, 205.13, 207.75, 201.78, 199.39, 200.58, 199.14, 192.45, 188.38, 191.24, 192.9, 192.18, 190.26, 187.88, 186.44, 183.09, 181.18, 194.75,
        194.75, 200.75, 203, 204, 207.5, 207.25, 209, 205.75, 207.5, 203.5, 202.25, 199.5, 202, 201, 198.25, 191.5, 187.75, 188.75, 190, 193, 193.75, 200, 200, 204, 194.5,
        192.25, 189.5, 188.5, 186.5, 187.5, 193.75, 196.5, 205, 204.25, 208, 211.75, 217, 213.25, 212.5, 213.75, 211.25, 214, 220.5, 212.5, 228, 225.75, 226.5, 233, 228.25, 225.5, 229, 229.5, 220.25, 220, 220, 223, 221, 216.5,
        211.25, 199.75, 192.5, 193.94, 192.74, 193.7, 195.6, 194.4, 195.36, 193.92, 185.79, 179.81, 162.13, 165.23, 168.81, 172.63, 168.81, 164.51, 159.01, 159.25, 154.95,
        152.08, 153.75, 153.74, 157.32, 160.66, 157.79, 152.54, 152.54, 156.84, 159.22, 162.08, 161.13, 166.38, 168.75, 73.75, 174.93, 182.8, 185.9, 213.58, 212.15,
        207.84, 214.76, 211.17, 206.62, 199.46, 196.37, 191.59, 183.95, 185.15, 172.97, 169.38, 160.08, 162.47, 159.6, 157.21, 147.91, 145.29, 146.01, 141.95, 143.38, 136, 124.55, 119.07, 116.69, 122.17, 122.41, 122.41, 136, 136.5,
        138.75, 138.5, 136.5, 135, 131.5, 130.75, 131.5, 133.25, 133.75, 134, 136.5, 133.75, 131, 121.5, 120.5, 115.5, 114.25, 110.25, 110, 110, 110.5, 111, 111.5, 110,
        112, 112, 113.25, 115.5, 120, 115.75, 112.75, 114, 112.25, 114.5, 117, 118.25, 127, 139.5, 131, 131.75, 133, 130, 127, 132.25, 128.5, 132.25, 134.5, 136.5, 135, 138.25,
        138, 134.5, 138, 138.5, 135, 130.75, 132.5, 129.75, 125.5, 124.5, 125, 120, 119, 121.05, 122.49, 122.97, 124.88, 119.61, 119.61, 126.79, 126.78, 127.26, 127.73, 127.26, 122.47, 118.17, 118.41, 118.65, 120.08, 123.42, 119.84,
        118.64, 122.71, 124.86, 121.27, 120.79
};




/*
 *  Reduction to produce the expected exposure profile EE[t] curve for the IRS
 * Expected Exposure  EPE(t) = 𝔼 [max(V , 0)|Ft]
*/
void reduce(std::vector<double>& expected_exposure, std::vector<std::vector<double>>& exposures, int timepoints_size, int simN) {
    for (int t = 0; t < timepoints_size; t++) {
        double sum = 0.0;
        for (int i = 0; i < simN; i++) {
            sum += std::max(exposures[i][t], 0.0);
        }
        expected_exposure[t] = sum /(double)simN;
    }
}

/*
 * Risk Factor Simulation Grid Matrix MAX_SIM x MAX_TENOR double array
 */
const int MAX_SIM = 10000;
const int MAX_TENOR = 51;
const double MAX_EXPIRY = 10.0;
const double DTAU = 0.5;

// Generate a grid of tenor_size x simN contains the risk factors Forward Rates F(t, T1, T2) from 0 .. T (expiry)
std::vector<std::vector<double>> fwd_rates(MAX_SIM, std::vector<double>(MAX_TENOR, 0.0));

// Generate one Exposure Profile by Simulation Step
std::vector<std::vector<double>> exposures(MAX_SIM, std::vector<double>(MAX_EXPIRY/DTAU, 0.0));

/*
 * Main Entry Point
// Initialize input parameters and trigger the mc_simulation to generate the averaged EE[t] for the IRS
// CVA is calculated using parameters PdF(t), DF(t) & EE(t)
// Expected exposure applied to a portfolio with only one trade IRS EURIBOR 6M 10Y
// Risk Factor Forward Rates F(t; t1, t2) is generated with HJM mdel
 */

int main() {
    const int tenors_size = 51;
    const double expiry = MAX_EXPIRY;
    const double dtau = 0.5;
    const int dimension = 3;

    // Calibration volatitliy & drift calibration using linear least square curve fitting
    //callibrate_volatilities_hjm(volatilities, drifts, yield_curve_EURIBOR_6M_10Y);

    // Discounts  or Zero Coupon Bond bootstrapped from the Spot Rates
    SpotRateYieldCurveTermStructure yieldCurve(spot_rates, 25.0, dtau);

    // Recovery Rates
    double recovery = 0.04;

    // Survival Probability Bootstrapping
    SurvivalProbabilityTermStructure survivalProbabilityCurve(tenors, spreads, yieldCurve, recovery, tenor_size); // dtau

    // Total Number of simulation points
    int simulation_points = expiry/dtau;

    // Expected Exposure Profile
    std::vector<double> expected_exposure(simulation_points, 0.0);

    // copy the spot rates as the first curve inside fwd_rates
    std::copy(spot_rates.begin(), spot_rates.end(), fwd_rates[0].begin());

    // Interest Rate Swap Instrument
    InterestRateSwap payOff(pricing_points, floating_schedule,  fixed_schedule, 10, 0.04, 10.0, 0.5);

    // Gaussian Random Number Generator
    ErfInvGaussianRandomGenerator erfinv;
    NormalRandomGenerator gaussian;

    // Simulation header ouput
    std::cout <<  "CVA" << "    " << "Iterations" << "    " << "Execution Time(s)" << std::endl;

    //  Increase the simulations numbers and analyze convergence and std deviation
    for (int simN = 500; simN < MAX_SIM; simN += 250) {
        // simulation step size
        double dt = expiry/simN;
        double duration = 0.0;

        // HJM Stochastic SDE Simulation Model
        HJMStochasticProcess hjm_sde(spot_rates, drifts, volatilities, fwd_rates, dimension, dt, dtau, tenors_size, 25.0);

        // Monte Carlo Simulation Engine generate the Exposure IRS Grid
        MonteCarloSimulation<ErfInvGaussianRandomGenerator, HJMStochasticProcess, InterestRateSwap> mc_engine(payOff, erfinv, hjm_sde, simN);
        mc_engine.calculate(exposures, duration);

        // Calculate Statistics max, median, quartiles, 97.5th percentile on exposures
        // Calculate Potential Future Exposure (PFE) at the 97.5th percentile and media of the Positive EE

        // Calculate Expected Exposure (EE)  EPE(t) = 𝔼 [max(V , 0)|Ft]
        // Use reduction to generate across timepoints the expected exposure profile for the IRS
        reduce(expected_exposure, exposures, simulation_points, simN);

        // Report Expected Exposure Profile Curve
        display_curve(expected_exposure);

        // Calculate The Unilateral CVA - Credit Value Adjustment Metrics Calculation.
        // For two conterparties A - B. Counterparty A want to know how much is loss in a contract due to B defaulting
        ExpectedExposureTermStructure expectedExposureCurve(pricing_points,expected_exposure, expiry);
        double cva = calculate_cva(recovery, yieldCurve, expectedExposureCurve, survivalProbabilityCurve, pricing_points, expiry);

        std::cout << std::setprecision(6)<< std::fixed << cva << " " << simN << " " << duration << std::endl;

        // Compare the exposures at maximum with percentiles consider
        // Sensitivity Analysis considering the very small and negative rates
    }

    exit(0);
}