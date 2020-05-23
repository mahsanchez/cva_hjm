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

using namespace std;

vector<double> tenors = {
        0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5,
        17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.5, 25.0 };

std::vector<double> timepoints = {
        0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5,
        17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.5, 25.0
};

/* Exposure Points
   daily calculations for 1 week
   weekly calculations up to 1 month
   biweekly up to 3 months
   monthly up to 1 year
   quarterly up to 5 years ]
   yearly up to the end time point up to 10Y
   Y5y Timesteps until 50Y
 */
/*  Convert days relatives to 30/360 calendar
    std::transform(exposure_points.begin(), exposure_points.end(), exposure_points.begin(), [](double day) {
        return day/360;
    });
*/
std::vector<double> exposure_timepoints = {
        1,2,3,4,5, 15, 22, 29, 30, 36, 50, 64, 76, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 780, 900, 1020, 1140, 1260, 1380, 1500, 1720, 1840, 2140, 2520, 2880, 3240, 3600, 4800, 7200
};

// Interest Rate Swap product 6M Euribor 10Y
std::vector<double> floating_schedule(timepoints);
std::vector<double> fixed_schedule = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0 };


void display_curve(std::vector<double> &curve) {
    std::cout << std::setprecision(6)<< std::fixed;
    std::copy(curve.begin(), curve.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
}

// TODO - Review Discount and Forward Rates
class YieldCurveTermStructure {
public:
    YieldCurveTermStructure(std::vector<double>& rates, std::vector<double>& tenors, int maturity) {
        interpolator.initialize(tenors, rates);
    }

    // discount factor at time t2 as seem from time t1
    double discount(double t1, double t2, double dt = 0.5) {
        double r = 0.0;
        for (double x = t1; x <= t2; x+= dt) {
            r += interpolator.find(x);
        }
        double result = std::exp(-r * dt);
        return result;
    }

    // Compute Forward Rate as L(t; 0.5, 1) = 1/CplTau*(EXP(MC!Caplet!C62*CplTau)-1)
    // L=m(e^(f/m)- 1), where m is compounding frequency per year  m = 1/0.5
    double forward_rate(double t1, double t2, double dtau) {
        double r = 1.0/dtau*(std::exp( interpolator.find(t1) *dtau)-1.0);
        return r;
    }

    /*
    double forward_rate(double t1, double t2, double dtau) {
        double r = 0.0;
        for (double t = t1; t <= t2; t += dtau) {
            r += points[t];
        }
        double result =  r / (t2/dtau + 1.0) ;
        return result;
    }
     */

private:
    linear_interpolator interpolator;
};


/*
 * VainillaInterestRateSwap
 */

class VanillaInterestRateSwap {
public:
    VanillaInterestRateSwap(double _notional, int initial_cashflow_, double t_, double _K, std::vector<double>& _floating_schedule, std::vector<double>& _fixed_schedule,  YieldCurveTermStructure & _yieldCurve) :
    notional(_notional), initial_cashflow(initial_cashflow_), K(_K), t(t_), floating_schedule(_floating_schedule), fixed_schedule(_fixed_schedule), yieldCurve(_yieldCurve)
    {
          calculate();
    };

    double floatingLeg() {
        double result = 0.0;
        for (int i = initial_cashflow; i < floating_schedule.size(); i++) {
            double t2 = floating_schedule[i];
            double t1 = floating_schedule[i-1];
            result += t2 * yieldCurve.forward_rate(t1, t2, t2 - t1) * yieldCurve.discount(t, t2, 0.5) ;
            //std::cout << "cash_flow " << i << " forward_rate " << yieldCurve.forward_rate(t1, t2, t2-t1) << " discount factor: " << yieldCurve.discount(t, t2, 0.5) << std::endl;
        }
        return result;
    }

    double fixedLeg() {
        double result = 0.0;
        for (int i = initial_cashflow; i < fixed_schedule.size(); i++) {
            double t2 = fixed_schedule[i];
            result += t2 * yieldCurve.discount(t, t2, 1.0) ;
        }
        result *= K;
        return result;
    }

    void calculate() {
        floating_leg = floatingLeg() ;
        fixed_leg = fixedLeg();
        npv = floating_leg - fixed_leg;
        npv *= notional;
    }

    double price() {
        return npv;
    }

private:
    int initial_cashflow;
    double notional;
    double t;
    double K;
    double floating_leg;
    double fixed_leg;
    double npv = 0.0;
    std::vector<double>& floating_schedule;
    std::vector<double>& fixed_schedule;
    YieldCurveTermStructure &yieldCurve;
};

/*
 * IRS Mark to Marking
*/
class InterestRateSwapExposureEngine {
public:
    InterestRateSwapExposureEngine(std::vector<double>& _exposures, std::vector<std::vector<double>>& _forward_rates, std::vector<double> timepoints_, double notional_, double K_, double maturity_, double dtau_, double dt_) :
    exposures(_exposures), forward_rates(_forward_rates), timepoints(timepoints_), notional(notional_), K(K_), maturity(maturity_), dtau(dtau_), dt(dt_){
        calculate();
    }

    void calculate() {
        int maturity = 51;
        for (int cashflow = 1; cashflow < maturity; cashflow++) {
            int t = timepoints[cashflow]/dt;
            YieldCurveTermStructure yieldCurve(forward_rates[t], timepoints, maturity);
            VanillaInterestRateSwap irs(notional, cashflow, timepoints[cashflow], K, floating_schedule, fixed_schedule, yieldCurve);
            exposures[cashflow] = irs.price();
        }
    }

private:
    double notional; // notional
    double K; // fixed rates IRS
    double maturity;
    double dtau;
    double dt;
    std::vector<double> timepoints;
    std::vector<double>& exposures;
    std::vector<std::vector<double>>& forward_rates;
};

/*
 * SurvivalProbabilityTermStructure
 * CDS Bootstrapping JPMorgan Methodology
 * VB code at http://mikejuniperhill.blogspot.com/2014/08/bootstrapping-default-probabilities.html
*/

class SurvivalProbabilityTermStructure {
public:
    SurvivalProbabilityTermStructure(std::vector<double>& timepoints, std::vector<double>& spreads, YieldCurveTermStructure& yieldCurve, double recovery, int maturity) : interpolator()
    {
        std::vector<double> probabilities(timepoints.size(), 0.0);
        bootstrap(probabilities, timepoints, spreads, yieldCurve, recovery, maturity);
        interpolator.initialize(timepoints, probabilities);
    }

    double operator() (double timepoint) {
        return interpolator.find(timepoint);
    }

    void bootstrap(std::vector<double>& p, std::vector<double>& timepoints, std::vector<double>& spreads, YieldCurveTermStructure& yieldCurve, double recovery, int maturity) {
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
                    term *= yieldCurve.discount(0.0, timepoints[j]);
                    terms += term;
                }

                double dtau = timepoints[i] - timepoints[i-1];
                divider = loss + dtau * spreads[i];
                divider *= yieldCurve.discount(0.0, timepoints[i]);
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
double calculate_cva(double recovery, YieldCurveTermStructure &yieldCurve, ExpectedExposureTermStructure &expected_exposure, SurvivalProbabilityTermStructure &survprob, std::vector<double> &exposure_points, int maturity, double dtau = 0.5) {
    double cva = 0.0;

    // TODO - USe parallel reduction to sum the vector array
    for (int i = 1; i < exposure_points.size(); i++) {
        double t = exposure_points[i];
        double t0 = exposure_points[i-1];
        cva += yieldCurve.discount(0.0, t) * expected_exposure(t) * (survprob(t0) - survprob(t) ) ;
    }

    cva = cva * (1 - recovery) ;

    return cva;
}



