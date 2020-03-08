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
        std::reverse(std::begin(rates), std::end(rates));

        std::vector<double> z(fwd_rates.size());
        std::partial_sum(rates.begin(), rates.end(), z.begin());
        std::reverse(std::begin(z), std::end(z));

        double t = 0;
        std::transform(z.begin(), z.end(), z.begin(), [&dtau](double r) {
            return std::exp( -(dtau * r) );
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
        for (int i = initial_cashflow + 1; i < floating_schedule.size(); i++) {
            double point = floating_schedule[i];
            result += dtau * L(point) * notional * dcf(point) ;
        }
        return result;
    }

    double fixedLeg() {
        double result = 0.0;
        double dtau = 0.5;
        for (int i = initial_cashflow + 1; i < fixed_schedule.size(); i++) {
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
    CVAPricer(double LGD, TermStructure eexposure, TermStructure pd, DiscountFactorsTermStructure dcf) :
            lgd(LGD), ee_curve(eexposure), pd_curve(pd), df_curve(dcf)
    {
        calculate();
    }

    inline double cva(double t) {
        double result = ee_curve(t);
        result *= pd_curve(t);
        result *= df_curve(t);
        return result;
    }

    // Apply Trapezoidal Rule for Integration
    void calculate() {
        double dtau = 0.5;
        // Computing sum of first and last terms in above formula
        double s = cva(0.0) + cva(25.0);
        // Adding middle terms in above formula
        for (int i = 1; i < 51; i++) {
            double delta = i * dtau;
            s += 2 * cva(i*dtau);
        }
        // h/2 indicates (b-a)/2n. Multiplying h/2 with s.
        value = 0.5 * dtau * s;
        // calculate CVA
        value = lgd * 0.5 * value;
    }

    double price() {
        return value;
    }

private:
    double value;
    double lgd;
    TermStructure ee_curve;
    TermStructure pd_curve;
    DiscountFactorsTermStructure df_curve;
};


