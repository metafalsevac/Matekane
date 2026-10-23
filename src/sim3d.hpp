// 3D six-degree-of-freedom rigid-body model with quaternion attitude.
//
// Frames (see README):
//   inertial  : flat-earth NED  (x North = down-valley distance, y East, z Down; altitude h = -z)
//   body      : FRD (x forward/nose, y right wing, z down), origin at the CG
//   attitude  : unit quaternion q = [w, x, y, z] (scalar first, Hamilton), body -> NED
//
// State (13): y = [ pN(3) | vB(3) = (u, v, w) | q(4) | omegaB(3) = (p, q, r) ]
//   pN_dot     = R(q) vB
//   vB_dot     = F_B / m - omegaB x vB
//   q_dot      = 1/2 q (x) [0, omegaB]
//   omegaB_dot = I^{-1} (M_B - omegaB x I omegaB)
// F_B = R(q)^T [0, 0, m g] + [T, 0, 0] + aero ; M_B = aerodynamic moments about the CG.
//
// Aerodynamics: a linear stability-derivative generalisation of the 2D model.
//   alpha = atan2(w, u), beta = asin(v / V), qbar = 1/2 rho V^2,
//   phat = p b/(2V), qhat = q c/(2V), rhat = r b/(2V)
//   CL = CL_cmd + CLa alpha + CLq qhat + CLde de     (CL_cmd = the 2D swept CL, i.e. the lift
//                                                     coefficient of the high-lift configuration
//                                                     at zero body angle of attack)
//   CD = CD0 + K CL^2                                 (same polar as 2D)
//   CY = CYb beta + CYdr dr
//   Cl = Clb beta + Clp phat + Clr rhat + Clda da + Cldr dr
//   Cm = Cm0 + Cma alpha + Cmq qhat + Cmde de
//   Cn = Cnb beta + Cnp phat + Cnr rhat + Cnda da + Cndr dr
// Lift/drag/side force act in wind axes and are rotated into body axes; thrust acts along +x_B.
// With Cm0 = 0 the aircraft is trimmed at alpha = 0 where CL = CL_cmd and thrust is aligned with
// the velocity, i.e. exactly the 2D model's assumptions. Constraining the motion to the pitch plane
// (v = p = r = 0, no lateral controls) and making the pitch axis stiff therefore reproduces the
// 2D results (verified in tests).
#pragma once
#include <array>
#include <cmath>

#include "model.hpp"
#include "vecmath.hpp"

namespace mk {

struct Params3D {
    double mass = MASS;
    // inertia tensor about CG in body axes, slug*ft^2 (representative single-engine light aircraft,
    // Cessna-172 class); Ixz is the only product of inertia allowed by left/right symmetry
    double Ixx = 948.0, Iyy = 1346.0, Izz = 1967.0, Ixz = 0.0;
    double S = S_REF;
    double b = 36.0;            // wing span, ft
    double cbar = S_REF / 36.0; // mean aerodynamic chord, ft (= S/b = 5.4375)
    double CL_cmd = 1.5;        // see header comment
    double CLa = 4.6, CLq = 3.9, CLde = 0.43;
    double CYb = -0.31, CYdr = 0.187;
    double Clb = -0.089, Clp = -0.47, Clr = 0.096, Clda = 0.178, Cldr = 0.0147;
    double Cm0 = 0.0, Cma = -0.89, Cmq = -12.4, Cmde = -1.28;
    double Cnb = 0.065, Cnp = -0.03, Cnr = -0.099, Cnda = -0.053, Cndr = -0.0657;
    double CL_max = CL_MAX;     // reported (stall is not modelled; see README)
    // control deflections, rad (held constant). Positive da -> right roll, positive de -> nose down,
    // positive dr -> nose left (standard sign conventions)
    double da = 0.0, de = 0.0, dr = 0.0;
    // switches (used by the sanity tests)
    bool aero = true, propulsion = true, gravity = true;
    // renormalise q after each non-Lie step (Sys3D::project)
    bool renormalize = true;
};

struct AeroState {
    double V, alpha, beta, CL, qbar;
};

struct Sys3D {
    static constexpr std::size_t N = 13;
    static constexpr bool has_lie = true;
    static constexpr std::size_t NL = 12;  // local chart: [pN | vB | theta(3) | omegaB]
    Params3D P;
    std::array<V3, 3> I{}, Iinv{};

    explicit Sys3D(const Params3D& p) : P(p) {
        I = {{{P.Ixx, 0.0, -P.Ixz}, {0.0, P.Iyy, 0.0}, {-P.Ixz, 0.0, P.Izz}}};
        // inverse of the symmetric (x-z coupled) inertia tensor
        const double det = P.Ixx * P.Izz - P.Ixz * P.Ixz;
        Iinv = {{{P.Izz / det, 0.0, P.Ixz / det}, {0.0, 1.0 / P.Iyy, 0.0}, {P.Ixz / det, 0.0, P.Ixx / det}}};
    }

    // ---- state packing helpers ----
    static V3 pos(const Vec<N>& y) { return {y[0], y[1], y[2]}; }
    static V3 vel(const Vec<N>& y) { return {y[3], y[4], y[5]}; }
    static Quat quat(const Vec<N>& y) { return {y[6], y[7], y[8], y[9]}; }
    static V3 omega(const Vec<N>& y) { return {y[10], y[11], y[12]}; }
    static Vec<N> pack(const V3& p, const V3& v, const Quat& q, const V3& w) {
        return {p[0], p[1], p[2], v[0], v[1], v[2], q.w, q.x, q.y, q.z, w[0], w[1], w[2]};
    }
    static double dist(const Vec<N>& y) { return y[0]; }
    static double alt(const Vec<N>& y) { return -y[2]; }
    static double airspeed(const Vec<N>& y) { return norm(vel(y)); }

    // Initial condition matching the 2D scenario: level flight path (gamma0) at V0, wings level,
    // heading north, alpha0 = 0 so theta0 = gamma0; optional initial pitch rate.
    static Vec<N> initial(double pitch_rate0 = 0.0, double alpha0 = 0.0) {
        Quat q = qfrom_euler(0.0, GAMMA0 + alpha0, 0.0);
        return pack({X0, 0.0, -H0}, {V0 * std::cos(alpha0), 0.0, V0 * std::sin(alpha0)}, q, {0.0, pitch_rate0, 0.0});
    }

    AeroState aero_state(const Vec<N>& y) const {
        const V3 v = vel(y);
        const double V = norm(v);
        AeroState a{V, 0.0, 0.0, P.CL_cmd, 0.5 * rho(alt(y)) * V * V};
        if (V > 1e-9) {
            a.alpha = std::atan2(v[2], v[0]);
            a.beta = std::asin(std::max(-1.0, std::min(1.0, v[1] / V)));
            const double qhat = y[11] * P.cbar / (2 * V);
            a.CL = P.CL_cmd + P.CLa * a.alpha + P.CLq * qhat + P.CLde * P.de;
        }
        return a;
    }

    // Forces and moments in body axes.
    void forces_moments(const Vec<N>& y, V3& F, V3& M) const {
        const Quat q = quat(y);
        const V3 w = omega(y);
        F = {0.0, 0.0, 0.0};
        M = {0.0, 0.0, 0.0};
        if (P.gravity) F = ned_to_body(q, {0.0, 0.0, P.mass * G});
        const AeroState a = aero_state(y);
        if (P.propulsion) F[0] += thrust(a.V, alt(y));
        if (P.aero && a.V > 1e-9) {
            const double V = a.V, al = a.alpha, be = a.beta;
            const double ph = w[0] * P.b / (2 * V), qh = w[1] * P.cbar / (2 * V), rh = w[2] * P.b / (2 * V);
            const double CL = a.CL;
            const double CD = CD0 + K_IND * CL * CL;
            const double CY = P.CYb * be + P.CYdr * P.dr;
            const double qS = a.qbar * P.S;
            const double D = qS * CD, Y = qS * CY, L = qS * CL;
            // wind -> body: F_B = C_bw [-D, Y, -L]
            const double ca = std::cos(al), sa = std::sin(al), cb = std::cos(be), sb = std::sin(be);
            F[0] += ca * cb * (-D) - ca * sb * Y - sa * (-L);
            F[1] += sb * (-D) + cb * Y;
            F[2] += sa * cb * (-D) - sa * sb * Y + ca * (-L);
            const double Cl = P.Clb * be + P.Clp * ph + P.Clr * rh + P.Clda * P.da + P.Cldr * P.dr;
            const double Cm = P.Cm0 + P.Cma * al + P.Cmq * qh + P.Cmde * P.de;
            const double Cn = P.Cnb * be + P.Cnp * ph + P.Cnr * rh + P.Cnda * P.da + P.Cndr * P.dr;
            M = {qS * P.b * Cl, qS * P.cbar * Cm, qS * P.b * Cn};
        }
    }

    // Translational + rotational rates excluding the quaternion part.
    void rates(const Vec<N>& y, V3& pdot, V3& vdot, V3& wdot) const {
        V3 F, M;
        forces_moments(y, F, M);
        const V3 v = vel(y), w = omega(y);
        pdot = body_to_ned(quat(y), v);
        const V3 wxv = cross(w, v);
        vdot = {F[0] / P.mass - wxv[0], F[1] / P.mass - wxv[1], F[2] / P.mass - wxv[2]};
        const V3 Iw = {dot(I[0], w), dot(I[1], w), dot(I[2], w)};
        const V3 wxIw = cross(w, Iw);
        const V3 rhs = {M[0] - wxIw[0], M[1] - wxIw[1], M[2] - wxIw[2]};
        wdot = {dot(Iinv[0], rhs), dot(Iinv[1], rhs), dot(Iinv[2], rhs)};
    }

    Vec<N> deriv(double /*t*/, const Vec<N>& y) const {
        V3 pd, vd, wd;
        rates(y, pd, vd, wd);
        const Quat q = quat(y);
        const V3 w = omega(y);
        const Quat qd = qmul(q, Quat{0.0, w[0], w[1], w[2]});
        return {pd[0], pd[1], pd[2], vd[0], vd[1], vd[2], 0.5 * qd.w, 0.5 * qd.x, 0.5 * qd.y, 0.5 * qd.z,
                wd[0], wd[1], wd[2]};
    }

    // Projection onto the unit sphere after each classical (non-Lie) step.
    void project(Vec<N>& y) const {
        if (!P.renormalize) return;
        const double n = std::sqrt(y[6] * y[6] + y[7] * y[7] + y[8] * y[8] + y[9] * y[9]);
        for (int i = 6; i < 10; ++i) y[i] /= n;
    }

    // ---- Lie-group (RKMK) chart around a base state: q = q_base (x) exp(theta) ----
    Vec<NL> lie_zero(const Vec<N>& base) const {
        return {base[0], base[1], base[2], base[3], base[4], base[5], 0.0, 0.0, 0.0, base[10], base[11], base[12]};
    }
    Vec<N> lie_to_global(const Vec<N>& base, const Vec<NL>& z) const {
        const Quat q = qmul(quat(base), qexp({z[6], z[7], z[8]}));
        return pack({z[0], z[1], z[2]}, {z[3], z[4], z[5]}, q, {z[9], z[10], z[11]});
    }
    Vec<NL> lie_deriv(double /*t*/, const Vec<N>& base, const Vec<NL>& z) const {
        const Vec<N> y = lie_to_global(base, z);
        V3 pd, vd, wd;
        rates(y, pd, vd, wd);
        const V3 thd = dexpinv_right({z[6], z[7], z[8]}, omega(y));
        return {pd[0], pd[1], pd[2], vd[0], vd[1], vd[2], thd[0], thd[1], thd[2], wd[0], wd[1], wd[2]};
    }
};

}  // namespace mk
