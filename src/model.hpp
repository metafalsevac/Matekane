// Shared aircraft, atmosphere and terrain model (US customary units: ft, slug, lbf, s, deg Rankine).
#pragma once
#include <cmath>

namespace mk {

// ---- atmosphere (ISA troposphere) ----
constexpr double T0 = 518.67;     // sea level temperature, R
constexpr double P0 = 2116.224;   // sea level pressure, lbf/ft^2
constexpr double LAPSE = -0.003566;  // temperature lapse rate, R/ft
constexpr double RGAS = 1716.0;   // gas constant of air, ft*lbf/(slug*R)
constexpr double P_EXP = 5.2561;  // ISA pressure exponent g0/(|lapse| R) (standard value)

// ---- aircraft (values from the original assignment) ----
constexpr double G = 32.2;        // gravitational acceleration, ft/s^2
constexpr double WEIGHT = 2712.5; // lbf (3100 lbf reduced by 12.5%)
constexpr double MASS = WEIGHT / G;  // slug (= 84.2391...)
constexpr double S_REF = 195.75;  // wing reference area, ft^2 (174 ft^2 + 12.5%)
constexpr double CL_MAX = 2.195;  // stall lift coefficient
constexpr double CD0 = 0.025;     // drag polar CD = CD0 + K_IND * CL^2
constexpr double K_IND = 0.056;
constexpr double T0SL = 714.9375; // max sea level static thrust, lbf (635.5 + 12.5%)
constexpr double A_SL = 0.0334;   // thrust airspeed coefficient, lbf*s^2/ft^2

// ---- scenario: takeoff from the Matekane airstrip ----
constexpr double V0 = 83.0;       // initial true airspeed, ft/s
constexpr double GAMMA0 = 0.0;    // initial flight path angle, rad
constexpr double X0 = 0.0;        // initial horizontal distance, ft
constexpr double H0 = 7649.63;    // initial altitude, ft (= Ground(0))
constexpr double X_LEVEL = 6233.58;  // ground is level beyond this distance
constexpr double H_LEVEL = 7480.58;  // level-ground elevation, ft
constexpr double X_END = 10000.0;    // stop when this distance is passed

inline double temp(double h) { return T0 + LAPSE * h; }
inline double press(double h) { return P0 * std::pow(1.0 + LAPSE * h / T0, P_EXP); }
inline double rho(double h) { return press(h) / (RGAS * temp(h)); }
inline double rho0() { return rho(0.0); }
inline double veas(double vtrue, double h) { return std::sqrt(rho(h) / rho0()) * vtrue; }

// Available thrust T = sigma * (T0SL - A_SL V^2), sigma = rho/rho0.
inline double thrust(double V, double h) { return rho(h) / rho0() * (T0SL - A_SL * V * V); }
inline double drag_coeff(double CL) { return CD0 + K_IND * CL * CL; }

// Terrain elevation (ft) along the valley as a function of horizontal distance x (ft).
inline double ground(double x) {
    if (x > X_LEVEL) return H_LEVEL;  // polynomial is continuous here (7480.579)
    return -7.54673e-12 * x * x * x * x + 9.16044e-8 * x * x * x - 1.8723e-4 * x * x - 5.91546e-1 * x + 7.64963e3;
}

// Number of values in the sweep min, min+step, ..., <= max. Uses an integer count because
// accumulating `x += step` in floating point can skip the last value.
inline int sweep_count(double min, double max, double step) {
    if (!(step > 0.0) || max < min) return 0;
    return static_cast<int>(std::floor((max - min) / step + 1e-9)) + 1;
}

}  // namespace mk
