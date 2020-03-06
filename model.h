#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <chrono>
#include <ctime>
#include <random>

using namespace std;

class StochasticProcess {
public:
    StochasticProcess(double dt_, int dimension_, double dtau_, vector<double>& drifts_, vector<double> &vols0, vector<double> &vols1, vector<double> &vols2) :
            dt(dt_) , dimension(dimension_), dtau(dtau_), drifts(drifts_), _vols0(vols0), _vols1(vols1), _vols2(vols2)
    {}

    void evolve(vector<double>& s0, vector<double>& s, double phi0, double phi1, double phi2) {
        for (int i = 0; i < dimension; i++) {
            double Mu = drifts[i] * dt;

            double vol = _vols0[i] * phi0;;
            vol += _vols1[i] * phi1;
            vol += _vols2[i] * phi2;
            double difussion = std::sqrt(dt) * vol;

            int t0;
            int t1;
            if (i < dimension - 1) {
                t0 = i;
                t1 = i+1;
            }
            else {
                t0 = i-1;
                t1 = i;
            }
            double dF = s0[t1] - s0[t0];

            double adjustment = (dF/dtau)* dt;

            s[i] = s0[i] + Mu + difussion + adjustment;
        }
    }

private:
    double dt;
    int dimension;
    double dtau;
    vector<double>& drifts;
    vector<double>& _vols0;
    vector<double>& _vols1;
    vector<double>& _vols2;
};





 /*
//template <class RNG, class S = Statistics>
// McSimulation(bool antitheticVariate, bool controlVariate);
extern "C" void MC_HJM_Musiela_CPU(std::vector<std::vector<double>> &grid, std::vector<double> &spot_rates, std::vector<double> &drifts, std::vector<std::vector<double>> &volatility, std::vector<std::vector<double>> &phi, const int tenors, double dt,  double tau)
{
    Statistics statistics;
    bool antitheticVariate;
    bool controlVariate;

    // MC path Generation
    //path_generator(grid, spot_rates, drifts, volatility, phi, tenors, dt, tau);
}
*/


// Intensity Model
// Default Probabilities

// fetch the credit spread information [51] timepoints from datagrapple.com