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
        return points.find(timepoint)->second;
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
        return L(points.find(timepoint)->second);
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
        return points.find(timepoint)->second;
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
    InterestRateSwap(double _notional, float _initial_cashflow, int _maturity, double _K, vector<double>& _floating_schedule, vector<double>& _fixed_schedule, YieldCurveTermStructure &forward_rate, DiscountFactorsTermStructure &dct_factors) :
    notional(_notional), initial_cashflow(_initial_cashflow), maturity(_maturity), K(_K), floating_schedule(_floating_schedule), fixed_schedule(_fixed_schedule), L(forward_rate), dcf_curve(dct_factors)
    {
    };

    double price() {
        double _fixedLeg = 0.0, _floatingLeg = 0.0;
        double date = 0;

        //Find the first table entry whose value is >= caller's x value
        auto iter = std::upper_bound(fixed_schedule.begin(), fixed_schedule.end(), initial_cashflow);
        if (iter != fixed_schedule.end()) {
            _fixedLeg = std::accumulate(iter, fixed_schedule.end(), 0.0, [this](double val, double x) {
                return val + (x * K * dcf_curve(x));
            });
        }

        iter = std::upper_bound(floating_schedule.begin(), floating_schedule.end(), initial_cashflow);
        if (iter != floating_schedule.end()) {
            _floatingLeg = std::accumulate(iter, floating_schedule.end(), 0.0, [this](double val, double x) {
                return val + (x * L(x) * dcf_curve(x));
            });
        }

        double npv = notional * (_floatingLeg - _fixedLeg);

        return npv;
    }

private:
    int initial_cashflow;
    double notional;
    double maturity;
    double K;
    std::vector<double>& floating_schedule;
    std::vector<double>& fixed_schedule;
    YieldCurveTermStructure& L;
    DiscountFactorsTermStructure &dcf_curve;
};



