#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <chrono>
#include <ctime>
#include <random>
#include <memory>
#include <cmath>

//#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include "q_numerics.h"

using namespace std;

// TODO - Implement an ExpectedExposureTermStructure (grid of mark to market) and return the EE(t) curve with interpolation
// TODO - Implement an ProbabilityOfDefaultTermStructure (grid of mark to market) and return the EE(t) curve with interpolation

std::vector<double> timepoints = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5,
                                  17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.5, 25.0 };

std::vector<double> floating_schedule(timepoints);

std::vector<double> fixed_schedule = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0 };


class TermStructure {
public:
    TermStructure(std::vector<double>& curve, double t0 = 0, double h = 0.5) {
        double x = 0.0;
        for ( auto iter = curve.begin(); iter != curve.end(); iter++) {
            points.insert({x, *iter});
            x += h;
        }
    }

    double operator() (double timepoint) {
        return points[timepoint];
    }

private:

    std::map<double,double,std::less<double>> points;
};


class YieldCurveTermStructure {
public:
    YieldCurveTermStructure(std::vector<double>& curve, double t0 = 0, double h = 0.5) {
        double x = 0.0;
        for ( auto iter = curve.begin(); iter != curve.end(); iter++) {
            points.insert({x, *iter});
            x += h;
        }
    }

    double operator() (double timepoint) {
        //return L(points[timepoint]);
        return points[timepoint];
    }

private:
    inline double L(double r, double tau = 0.5) {
        return (1/tau) * ( std::exp(r * tau) - 1 );
    }

    std::map<double,double,std::less<double>> points;
};


class DiscountFactorsTermStructure {
public:
    DiscountFactorsTermStructure(std::vector<double>& fwd_rates, double dtau = 0.5)
    {
        std::vector<double> rates(fwd_rates);
        std::vector<double> z(fwd_rates.size());

        std::partial_sum(rates.begin(), rates.end(), z.begin());
        std::transform(z.begin(), z.end(), z.begin(), [&dtau](double r) {
            double discount = std::exp( -r * dtau );
            return discount;
        });

        double x = 0.0;
        for ( auto iter = z.begin(); iter != z.end(); iter++) {
            points.insert({x, *iter});
            x += dtau;
        }
    }

    double operator() (double timepoint) {
        return points[timepoint];
    }

private:
    std::map<double,double,std::less<double>> points;
};


class DefaultProbabilityTermStructure {
public:
    DefaultProbabilityTermStructure(std::vector<double> &timepoints, TermStructure &spreads, DiscountFactorsTermStructure &discountFactorCurve,  double L, double dtau) {
        const int N = timepoints.size();

        std:vector<double> psurv(N, 1.0);

        psurv[1] = L / (spreads(timepoints[1])*dtau + L);

        points[0.0] = 1.0;
        points[0.5] = psurv[1];

        for (int i = 2; i < N; i++) {
            double timepoint = timepoints[i];
            double spread = spreads(timepoint);
            double psurvival = 0.0;
            for (int t = 1; t < i; t++) {
                psurvival += L * psurv[t-1];
                psurvival -= (L + dtau * spread);
                psurvival *= psurv[t];
                psurvival *= discountFactorCurve(timepoints[t]);
            }
            psurvival /= ( discountFactorCurve(timepoint) * (L + dtau * spread) );
            psurvival += (psurv[i-1] * L) / ( L + dtau * spread);
            psurv[i] = psurvival;
            points[timepoint] = psurvival > 0 ? psurvival : -psurvival;
        }
    }

    double operator() (double timepoint) {
        return points[timepoint];
    }

private:

    std::map<double,double,std::less<double>> points;
};

/*
* Given MC HJM simulation Grid extract Forward Rate and ZCB
        Settlement Days
        DayCounter
ReferenceDay
        Calendar

DayCount 30/360 Convention  ( 360*(y2 -y1) + 30*(m2 - m1) + (d2 - d1) ) / 360

 USes a cubic spline interpolator to obtain the values
*/

// TODO - Compute Forward Rate as L(t; 0.5, 1) = 1/CplTau*(EXP(MC!Caplet!C62*CplTau)-1)
// TODO L=m(e^(f/m)- 1), where m is compounding frequency per year  m = 1/0.5

/*
 * Vainilla Interest Rate Swap valuation
nominal: 1M
maturity: 20 nov 2022
init date: 20 nov 2012
schedule: 3m euribor
doubleing leg pays: Euribor 3m + 2% every 3 months
fixed leg pays:
 4% anually
last fixing for Euribor 3m : 1%
calendar Target
day counting convention fixed leg: 30/360
day counting convention doubleing leg: Actual/360

Market Data Forward Rate Curve **
Market Data Discount Factor Curve **

const std::vector<Date>& dates,
const std::vector<double>& dates,
const std::vector<double>& yields,
*/
class InterestRateSwap {
public:
    InterestRateSwap(double _notional, int start_day, int _maturity, double _K, vector<double>& _floating_schedule, vector<double>& _fixed_schedule, YieldCurveTermStructure &forward_rate, DiscountFactorsTermStructure &dct_factors) :
    notional(_notional), initial_cashflow(start_day), maturity(_maturity), K(_K), floating_schedule(_floating_schedule), fixed_schedule(_fixed_schedule), L(forward_rate), dcf(dct_factors)
    {
          calculate();
    };

    double floatingLeg() {
        double result = 0.0;
        double dtau = 0.5;
        for (int i = initial_cashflow; i < floating_schedule.size(); i++) {
            double point = floating_schedule[i];
            result += dtau * L(point) * notional * dcf(point) ;
        }
        return result;
    }

    double fixedLeg() {
        double result = 0.0;
        double dtau = 0.5;
        for (int i = initial_cashflow; i < fixed_schedule.size(); i++) {
            double point = fixed_schedule[i];
            result += dtau * K * notional * dcf(point) ;
        }
        return result;
    }

    void calculate() {
        floating_leg = floatingLeg() ;
        fixed_leg = fixedLeg();
        npv = floating_leg - fixed_leg;
    }

    double price() {
        return npv;
    }

private:
    int initial_cashflow;
    double notional;
    double maturity;
    double K;
    double floating_leg;
    double fixed_leg;
    double npv = 0.0;
    std::vector<double>& floating_schedule;
    std::vector<double>& fixed_schedule;
    YieldCurveTermStructure& L;
    DiscountFactorsTermStructure &dcf;
};


// CVA pricing using trapezoidal method

class CVAPricer {
public:
    CVAPricer(double LGD, TermStructure eexposure, DefaultProbabilityTermStructure pd, DiscountFactorsTermStructure dcf) :
            cva(0.0),lgd(LGD), ee_curve(eexposure), pd_curve(pd), df_curve(dcf)
    {
        calculate();
    }

    // Apply Trapezoidal Rule for Integration
    void calculate() {
        double dtau = 0.5;
        for (int i = 1; i < 51; i++) {
            float t = timepoints[i];
            float t0 = timepoints[i-1];
            double value = ee_curve(t) - ee_curve(t0);
            value *= pd_curve(t0) - pd_curve(t);
            cva += value;
        }
        cva *= -lgd;
    }

    double price() {
        return cva;
    }

private:
    double cva;
    double lgd;
    TermStructure ee_curve;
    DefaultProbabilityTermStructure pd_curve;
    DiscountFactorsTermStructure df_curve;
};


