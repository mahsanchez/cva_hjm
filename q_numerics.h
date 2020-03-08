#include <utility>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <map>

// TODO - implement own integration/interpolation keep it simple using tra[ezium rule

template <typename F>
float trapezoidal(F f, double a, double b, int n, double h)
{
    // Computing sum of first and last terms in above formula
    double s = f(a) + f(b);

    // Adding middle terms in above formula
    for (int i = 1; i < n; i++)
        s += 2 * f(a + i*h);

    // h/2 indicates (b-a)/2n. Multiplying h/2 with s.
    return 0.5 * h * s;
}


class linear_interpolator {
public:
    linear_interpolator(const std::vector<double>& values, double t0 = 0.0, double h = 0.5)
    {
        double x = 0.0;
        double index = 0;

        for (int i = 0; i < t0; i++) {
            points.push_back( {x, 0.0} );
        }

        for (int i = t0; i < values.size(); i++) {
            x = index * h; index++;
            points.push_back( {x, values[i]} );
        }
    }

    double operator() (double x) const {
        if (points.size() == 0) {
            return 0.0;
        }

        //Define a lambda that returns true if the x value of a point pair is < the caller's x value
        auto lessThan =
                [](const std::pair<double, double>& point, double x)
                {return point.first < x;};

        //Find the first table entry whose value is >= caller's x value
        auto iter =
                std::lower_bound(points.cbegin(), points.cend(), x, lessThan);

        //If the caller's X value is greater than the largest
        //X value in the table, we can't interpolate.
        if(iter == points.cend()) {
            return (points.cend() - 1)->second;
        }

        //If the caller's X value is less than the smallest X value in the table,
        //we can't interpolate.
        if(iter == points.cbegin() and x <= points.cbegin()->first) {
            return points.cbegin()->second;
        }

        //We can interpolate!
        double upperX{iter->first};
        double upperY{iter->second};
        double lowerX{(iter - 1)->first};
        double lowerY{(iter - 1)->second};

        double deltaY{upperY - lowerY};
        double deltaX{upperX - lowerX};

        return lowerY + ((x - lowerX)/ deltaX) * deltaY;
    }

private:
    std::vector<std::pair<double, double>> points;
};
