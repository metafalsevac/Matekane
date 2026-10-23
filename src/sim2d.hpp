// 2D (vertical plane) point-mass model in wind axes -- the original Matekane formulation.
//   state y = [V, gamma, h, x]  (true airspeed ft/s, flight-path angle rad, altitude ft, distance ft)
//   V_dot     = (T - D)/m - g sin(gamma)            (thrust assumed along the velocity vector)
//   gamma_dot = L/(m V) - g cos(gamma)/V
//   h_dot     = V sin(gamma)
//   x_dot     = V cos(gamma)
// with L = 1/2 rho V^2 S CL, D = 1/2 rho V^2 S CD, CL held constant (the swept parameter).
#pragma once
#include "model.hpp"
#include "vecmath.hpp"

namespace mk {

struct Sys2D {
    static constexpr std::size_t N = 4;
    static constexpr bool has_lie = false;
    static constexpr std::size_t NL = 1;  // unused
    double CL = 1.5;
    double CD = drag_coeff(1.5);

    explicit Sys2D(double cl) : CL(cl), CD(drag_coeff(cl)) {}

    Vec<N> deriv(double /*t*/, const Vec<N>& y) const {
        const double V = y[0], gam = y[1], h = y[2];
        const double r = rho(h);
        const double Vdot = thrust(V, h) / MASS - r * V * V * S_REF * CD / (2 * MASS) - G * std::sin(gam);
        const double gdot = r * V * S_REF * CL / (2 * MASS) - G * std::cos(gam) / V;
        return {Vdot, gdot, V * std::sin(gam), V * std::cos(gam)};
    }
    void project(Vec<N>&) const {}
    // unused Lie hooks (never called; has_lie = false)
    Vec<NL> lie_zero(const Vec<N>&) const { return {0.0}; }
    Vec<N> lie_to_global(const Vec<N>& b, const Vec<NL>&) const { return b; }
    Vec<NL> lie_deriv(double, const Vec<N>&, const Vec<NL>&) const { return {0.0}; }

    static Vec<N> initial() { return {V0, GAMMA0, H0, X0}; }
    static double dist(const Vec<N>& y) { return y[3]; }
    static double alt(const Vec<N>& y) { return y[2]; }
    static double airspeed(const Vec<N>& y) { return y[0]; }
};

}  // namespace mk
