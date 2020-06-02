#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <chrono>
#include <ctime>
#include <memory>
#include <cmath>
#include <iostream>

#include "q_numerics.h"

#define DEBUG 100

using namespace std;

// Tenors as seen in the Spot Rate Curve coming from Bank of England
vector<double> tenors = {
        0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5,
        17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.5, 25.0 };

std::vector<double> timepoints = {
        0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5,
        17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.5, 25.0
};

// Cache Flow Schedule Interest Rate Swap product 6M Euribor 10Y
std::vector<double> pricing_points = {
        0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0
};

// Interest Rate Swap Cash Flows

std::vector<double> floating_schedule = {
        0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0
};

std::vector<double> fixed_schedule = {
        0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0
};

// Exposure Points
std::vector<double> exposure_timepoints = {
        1,2,3,4,5, 15, 22, 29, 30, 36, 50, 64, 76, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 780, 900, 1020, 1140, 1260, 1380, 1500, 1720, 1840, 2140, 2520, 2880, 3240, 3600, 4800, 7200
};


void display_curve(std::vector<double> &curve) {
    std::cout << std::setprecision(10)<< std::fixed;
    std::copy(curve.begin(), curve.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
}

void display_curve(std::vector<double> &curve, int begin, int end) {
    std::cout << std::setprecision(10)<< std::fixed;
    std::copy(curve.begin() + begin, curve.begin() + end, std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
}

void display_curve(std::vector<std::pair<double, double>> &curve) {
    std::cout << std::setprecision(10)<< std::fixed;
    for (int i = 0; i < curve.size(); i++) {
        std::pair<double, double> point = curve[i];
        std::cout << point.first << " " << point.second << std::endl;
    }
}


/*
 * Pricing Instrument Interest Rate Swap IRS
 */
struct InterestRateSwap {
    InterestRateSwap(std::vector<double> &pricing_points_, std::vector<double> &floating_schedule_,  std::vector<double> &fixed_schedule_, double notional_, double K_, double expiry_, double dtau_) :
            pricing_points(pricing_points_), floating_schedule(floating_schedule_), fixed_schedule(fixed_schedule_), notional(notional_), K(K_), expiry(expiry_), dtau(dtau_)
    {}

    std::vector<double> &pricing_points;
    std::vector<double> &floating_schedule;
    std::vector<double> &fixed_schedule;
    double notional;
    double K;
    double dtau;
    double expiry;
};


/*
 * The relationship between the discount function and the annually compounded yield curve, using a day count convention that reflects the
 * actual time between time t0 and t measured in years, can be written as
 * Reference day for a discount factor over the spot rate is today
 */
class SpotRateYieldCurveTermStructure {
public:
    SpotRateYieldCurveTermStructure(std::vector<double>& rates, double expiry, double dtau = 0.5) {
        // Build Discount Curve
        bootstrapDiscountFactorsCurve(rates, expiry, dtau);
    }

    // discount factor at time t2 as seem from time t1
    double discount(double t) {
        double df = discountCurve.find(t);
        return df;
    }

    void bootstrapDiscountFactorsCurve(std::vector<double>& rates, double expiry, double dtau = 0.5) {
        double tenor_size = expiry/dtau + 1;
        double tenor = 0.5;

        std::vector<double> partial_rates(rates.size(), 0.0);
        prefix_sum(&partial_rates[0], &rates[0], rates.size());

        discountCurve.add(0.0, 1.0);
        for (int i = 0; i < rates.size(); i++) {
            double discount = std::exp( -partial_rates[i] * dtau);
            discountCurve.add(tenor, discount);
            tenor += dtau;
        }
    }

private:
    linear_interpolator discountCurve;
};


class HJMYieldCurveTermStructure {
public:
    HJMYieldCurveTermStructure(std::vector<std::vector<double>>& fwds, double reference_day, double expiry, double dt_, double dtau_ = 0.5) {
        dtau = dtau_;
        dt = dt_;
        // Build Forward Rate
        bootstrapForwardCurve(fwds, reference_day, expiry, dt, dtau);
        // Build Discount Curve
        bootstrapDiscountFactorsCurve(fwds, reference_day, expiry, dt, dtau);
    }

    // discount factor at time t2 as seem from time t1
    double discount(double t) {
        double df = discountCurve.find(t);
        return df;
    }

    // Forward Libor Rate as  factor at time t2 as seem L(t; t1, t2)
    double forward(double t, double t2) {
        double df = forwardCurve.find(t);
        return df;
    }

    // Compute Forward LIBOR as L(t; t, T) At Fixed Tenors L=m(e^(f/m)- 1), where m is compounding frequency per year  m = 1/0.5
    void bootstrapForwardCurve(std::vector<std::vector<double>>& fwdGrid, double reference_day, double expiry, double dt, double dtau) {
        int simulation = reference_day/dt;
        int tenor_size = expiry/dtau;
        double tenor = 0.0;

        for (int i = 0; i < tenor_size; ++i) {
            double r = fwdGrid[simulation][i];
            double libor = 1.0/dtau*(std::exp(r*dtau)-1.0);
            forwardCurve.add(tenor, libor);
            tenor += dtau;
        }
    }

    /*
     * column to find the tenor related with the start day [ as seen from today [0], as seen in 6 months [0.5]  and so on
     * Entire row from fwdGrid begin to end; but only sum on tenors interested points tenor mod 50  dtau/dt
     */
    void bootstrapDiscountFactorsCurve(std::vector<std::vector<double>>& fwdGrid, double reference_day, double expiry, double dt, double dtau) {
        int tenor = reference_day/dtau; //
        int points_size = expiry/dtau;

        discountCurve.add(0.0, 1.0);

        int begin = 0;
        int leap = dtau/dt;
        int end = leap;

        double t = 0.0;
        for (int i = 1; i < points_size; i++) {
            double sum = 0.0;
            for (int sim = 0; sim < end; sim++) {
                sum += fwdGrid[sim][tenor];
            }
            double zcb = std::exp(-sum * dt);
            discountCurve.add(t, zcb);
            t += dtau;
            end += leap;
        }
        //Display Discount Curve
        //display_curve(discountCurve.getPoints());
    }

    linear_interpolator getDiscountCurve() {
        return discountCurve;
    }

private:

    double dtau;
    double dt;
    linear_interpolator forwardCurve;
    linear_interpolator discountCurve;
};


/*
 * VainillaInterestRateSwap
 * Valuation = Sum[0..T] (S - K) *  𝛼j * DF(t0, tj)
 *
 * 􏰷 𝛼j is the day-count fraction over which the coupon accrues.
􏰷 * B(t0 , tj ) is the value (at t0 ) of a discount factor maturing at time tj .
􏰷 * K is the fixed rate of the swap.
 * 􏰷S is the par swap rate where
 *
􏰷 The fixed and floating leg frequencies and day count bases are assumed to be the same for simplicity.
 */

class VanillaInterestRateSwapPricer { //InterestRateSwap
public:
    VanillaInterestRateSwapPricer(int initial_cashflow_, double reference_day_, InterestRateSwap& irs_,  HJMYieldCurveTermStructure & _yieldCurve) :
    initial_cashflow(initial_cashflow_), reference_day(reference_day_), irs(irs_), yieldCurve(_yieldCurve)
    {
          calculate();
    };

    void calculate() {
        double floatingL = floatingLeg();
        double fixedL = fixedLeg();
        npv = irs.notional * (floatingL - fixedL);
    }

    double floatingLeg() {
        double price = 0.0;
        for (int i = initial_cashflow; i < irs.floating_schedule.size(); i++) {
            double tau = 0.5;
            double t2 = floating_schedule[i];
            double t1 = floating_schedule[i-1];
            double sum = yieldCurve.forward(t1, t2);
            sum *= t2;
            sum *= yieldCurve.discount(t2);
            price += sum;
        }
       return price;
    }

    double fixedLeg() {
        double price = 0.0;
        for (int i = initial_cashflow; i < irs.fixed_schedule.size(); i++) {
            double t = fixed_schedule[i];
            double tau = 1.0;
            double sum = tau;
            sum *= irs.K;
            sum *= yieldCurve.discount(t);
            price += sum;
        }
        return price;
    }

    double price() {
        return npv;
    }

private:
    InterestRateSwap& irs;
    int initial_cashflow;
    double reference_day;
    double npv = 0.0;
    HJMYieldCurveTermStructure &yieldCurve;
};


/*
 * SurvivalProbabilityTermStructure
 * CDS Bootstrapping JPMorgan Methodology
 * VB code at http://mikejuniperhill.blogspot.com/2014/08/bootstrapping-default-probabilities.html
*/

class SurvivalProbabilityTermStructure {
public:
    SurvivalProbabilityTermStructure(std::vector<double>& timepoints, std::vector<double>& spreads, SpotRateYieldCurveTermStructure& yieldCurve, double recovery, int maturity) : interpolator()
    {
        std::vector<double> probabilities(timepoints.size(), 0.0);
        bootstrap(probabilities, timepoints, spreads, yieldCurve, recovery, maturity);
        interpolator.initialize(timepoints, probabilities);
    }

    double operator() (double timepoint) {
        return interpolator.find(timepoint);
    }

    void bootstrap(std::vector<double>& p, std::vector<double>& timepoints, std::vector<double>& spreads, SpotRateYieldCurveTermStructure& yieldCurve, double recovery, int maturity) {
        double loss = 1 - recovery;
        double term, terms, divider, term1, term2;

        std::transform(spreads.begin(), spreads.end(), spreads.begin(), [](double &s) {
            return s * 0.0001;
        });

        for (int i = 0; i < maturity; i++) {
            if (i == 0) {
                p[0] = 1.0;
            }
            else if (i == 1) {
                p[1] = loss / (spreads[1] * (timepoints[1] - timepoints[0]) + loss);
            }
            else {
                terms = 0.0;
                for (int j = 1; j < i; j++) {
                    double dtau = timepoints[j] - timepoints[j-1];
                    term = loss * p[j-1];
                    term -= (loss + dtau * spreads[i]) * p[j];
                    term *= yieldCurve.discount(timepoints[j]);
                    terms += term;
                }

                double dtau = timepoints[i] - timepoints[i-1];
                divider = loss + dtau * spreads[i];
                divider *= yieldCurve.discount(timepoints[i]);
                term1 = terms/divider;
                term2 = p[i-1] * loss;
                term2 /= (loss + dtau * spreads[i]);
                p[i] = term1 + term2;
            }
        }
    }
private:
    linear_interpolator interpolator;
};


/**
 * Expected Exposure Interpolated Curve
 */
class ExpectedExposureTermStructure {
public:
    ExpectedExposureTermStructure(std::vector<double>& timepoints, std::vector<double>& exposure_curve, int maturity) : interpolator() {
        interpolator.initialize(timepoints, exposure_curve);
    }

    double operator() (double timepoint) {
        return interpolator(timepoint);
    }

private:
    linear_interpolator interpolator;
};


/*
 * CVA Calculation
 * CVA =  E [ (1 - R) [ DF[t] * EE[t] * dPD[t] ] ]
 */
double calculate_cva(double recovery, SpotRateYieldCurveTermStructure &yieldCurve, ExpectedExposureTermStructure &expected_exposure, SurvivalProbabilityTermStructure &survprob, std::vector<double> &exposure_points, int maturity, double dtau = 0.5) {
    double cva = 0.0;

    for (int i = 1; i < exposure_points.size() - 1; i++) {
        double t = exposure_points[i];
        double t0 = exposure_points[i-1];
        cva += yieldCurve.discount(t) * expected_exposure(t) * (survprob(t0) - survprob(t) ) ;
    }

    cva = cva * (1 - recovery) ;

    return cva;
}



