// Time integrators for autonomous/non-autonomous ODE systems y' = f(t, y) on fixed-size vectors.
//
//   rk4        classic 4th-order Runge-Kutta, fixed step (default; same scheme as the original code)
//   rk45       Dormand-Prince 5(4) embedded pair, adaptive step size, FSAL, 4th-order dense output
//   lie-rk4    Runge-Kutta-Munthe-Kaas (RKMK) with the classic RK4 tableau
//   lie-rk45   RKMK with the Dormand-Prince tableau (adaptive)
//
// The Lie-group ("geometric") variants require the system to provide a local chart around the
// current state (see Sys3D): the attitude is advanced as q_{n+1} = q_n (x) exp(theta) where the
// rotation vector theta is integrated with the RK tableau from theta = 0 using
// theta_dot = dexp^{-1}(theta) * omega. Because exp() and quaternion products of unit quaternions
// are unit quaternions, |q| = 1 is preserved to rounding error with no renormalisation.
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "vecmath.hpp"

namespace mk {

enum class Method { RK4, RK45, LieRK4, LieRK45 };

inline Method parse_method(const std::string& s) {
    if (s == "rk4") return Method::RK4;
    if (s == "rk45" || s == "dopri5" || s == "dp45") return Method::RK45;
    if (s == "lie-rk4" || s == "rkmk4") return Method::LieRK4;
    if (s == "lie-rk45" || s == "rkmk45") return Method::LieRK45;
    throw std::invalid_argument("unknown integrator '" + s + "' (use rk4, rk45, lie-rk4, lie-rk45)");
}
inline const char* method_name(Method m) {
    switch (m) {
        case Method::RK4: return "rk4";
        case Method::RK45: return "rk45";
        case Method::LieRK4: return "lie-rk4";
        case Method::LieRK45: return "lie-rk45";
    }
    return "?";
}
inline bool is_lie(Method m) { return m == Method::LieRK4 || m == Method::LieRK45; }
inline bool is_adaptive(Method m) { return m == Method::RK45 || m == Method::LieRK45; }

// ----------------------------------------------------------------------------------------------
// Classic RK4: k1 = f(t, y); k2 = f(t+h/2, y+h/2 k1); k3 = f(t+h/2, y+h/2 k2); k4 = f(t+h, y+h k3)
//              y+ = y + h/6 (k1 + 2 k2 + 2 k3 + k4)
// ----------------------------------------------------------------------------------------------
template <std::size_t N, class F>
Vec<N> rk4_step(F&& f, double t, const Vec<N>& y, double h) {
    const Vec<N> k1 = f(t, y);
    const Vec<N> k2 = f(t + 0.5 * h, y + (0.5 * h) * k1);
    const Vec<N> k3 = f(t + 0.5 * h, y + (0.5 * h) * k2);
    const Vec<N> k4 = f(t + h, y + h * k3);
    Vec<N> r;
    for (std::size_t i = 0; i < N; ++i) r[i] = y[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    return r;
}

// ----------------------------------------------------------------------------------------------
// Dormand-Prince 5(4) (Hairer, Norsett & Wanner, "Solving ODEs I", DOPRI5).
// ----------------------------------------------------------------------------------------------
namespace dp {
constexpr double c2 = 1.0 / 5, c3 = 3.0 / 10, c4 = 4.0 / 5, c5 = 8.0 / 9;
constexpr double a21 = 1.0 / 5;
constexpr double a31 = 3.0 / 40, a32 = 9.0 / 40;
constexpr double a41 = 44.0 / 45, a42 = -56.0 / 15, a43 = 32.0 / 9;
constexpr double a51 = 19372.0 / 6561, a52 = -25360.0 / 2187, a53 = 64448.0 / 6561, a54 = -212.0 / 729;
constexpr double a61 = 9017.0 / 3168, a62 = -355.0 / 33, a63 = 46732.0 / 5247, a64 = 49.0 / 176,
                 a65 = -5103.0 / 18656;
constexpr double a71 = 35.0 / 384, a73 = 500.0 / 1113, a74 = 125.0 / 192, a75 = -2187.0 / 6784,
                 a76 = 11.0 / 84;  // = 5th-order weights b_i
// error coefficients e_i = b_i - bhat_i
constexpr double e1 = 71.0 / 57600, e3 = -71.0 / 16695, e4 = 71.0 / 1920, e5 = -17253.0 / 339200,
                 e6 = 22.0 / 525, e7 = -1.0 / 40;
// dense output (continuous extension of order 4)
constexpr double d1 = -12715105075.0 / 11282082432, d3 = 87487479700.0 / 32700410799,
                 d4 = -10690763975.0 / 1880347072, d5 = 701980252875.0 / 199316789632,
                 d6 = -1453857185.0 / 822651844, d7 = 69997945.0 / 29380423;
}  // namespace dp

template <std::size_t N>
struct DopriStep {
    Vec<N> y0, y1, k7;
    double err = 0.0;   // scaled RMS error norm (accept if <= 1)
    Vec<N> r2, r3, r4, r5;  // dense output coefficients
    // Continuous extension at t0 + theta*h, theta in [0, 1].
    Vec<N> dense(double theta) const {
        const double t1 = 1.0 - theta;
        Vec<N> r;
        for (std::size_t i = 0; i < N; ++i)
            r[i] = y0[i] + theta * (r2[i] + t1 * (r3[i] + theta * (r4[i] + t1 * r5[i])));
        return r;
    }
};

template <std::size_t N, class F>
DopriStep<N> dopri_step(F&& f, double t, const Vec<N>& y, const Vec<N>& k1, double h, double rtol,
                        double atol) {
    using namespace dp;
    Vec<N> tmp;
    for (std::size_t i = 0; i < N; ++i) tmp[i] = y[i] + h * a21 * k1[i];
    const Vec<N> k2 = f(t + c2 * h, tmp);
    for (std::size_t i = 0; i < N; ++i) tmp[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]);
    const Vec<N> k3 = f(t + c3 * h, tmp);
    for (std::size_t i = 0; i < N; ++i) tmp[i] = y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
    const Vec<N> k4 = f(t + c4 * h, tmp);
    for (std::size_t i = 0; i < N; ++i)
        tmp[i] = y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
    const Vec<N> k5 = f(t + c5 * h, tmp);
    for (std::size_t i = 0; i < N; ++i)
        tmp[i] = y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
    const Vec<N> k6 = f(t + h, tmp);
    DopriStep<N> s;
    s.y0 = y;
    for (std::size_t i = 0; i < N; ++i)
        s.y1[i] = y[i] + h * (a71 * k1[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
    s.k7 = f(t + h, s.y1);
    double acc = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        double e = h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * s.k7[i]);
        double sc = atol + rtol * std::max(std::fabs(y[i]), std::fabs(s.y1[i]));
        acc += (e / sc) * (e / sc);
    }
    s.err = std::sqrt(acc / N);
    for (std::size_t i = 0; i < N; ++i) {
        s.r2[i] = s.y1[i] - y[i];
        s.r3[i] = h * k1[i] - s.r2[i];
        s.r4[i] = s.r2[i] - h * s.k7[i] - s.r3[i];
        s.r5[i] = h * (d1 * k1[i] + d3 * k3[i] + d4 * k4[i] + d5 * k5[i] + d6 * k6[i] + d7 * s.k7[i]);
    }
    return s;
}

// ----------------------------------------------------------------------------------------------
// Generic driver.
// ----------------------------------------------------------------------------------------------
struct IntegratorOptions {
    double dt = 0.1;      // fixed step (rk4 / lie-rk4) and output interval (all methods)
    double rtol = 1e-8;   // adaptive methods
    double atol = 1e-8;
    double h0 = 0.0;      // initial step for adaptive methods (0 -> dt)
    double hmax = 0.0;    // max step for adaptive methods (0 -> unlimited)
    double hmin = 1e-12;
    long max_steps = 50000000;
    bool project = true;  // apply Sys::project (e.g. quaternion renormalisation) after non-Lie steps
};

struct IntegratorStats {
    long steps = 0;      // accepted steps
    long rejected = 0;   // rejected steps (adaptive)
    long fevals = 0;     // right-hand-side evaluations
};

// One fixed step of rk4 / lie-rk4.
template <class Sys>
Vec<Sys::N> fixed_step(const Sys& sys, Method m, double t, const Vec<Sys::N>& y, double h,
                       IntegratorStats& st, bool project) {
    constexpr std::size_t N = Sys::N;
    if (m == Method::LieRK4) {
        if constexpr (Sys::has_lie) {
            auto fl = [&](double tt, const Vec<Sys::NL>& z) { ++st.fevals; return sys.lie_deriv(tt, y, z); };
            Vec<Sys::NL> z1 = rk4_step(fl, t, sys.lie_zero(y), h);
            ++st.steps;
            return sys.lie_to_global(y, z1);
        }
        throw std::invalid_argument("Lie-group integrators need a system with a rotational state (3D mode)");
    }
    auto f = [&](double tt, const Vec<N>& yy) { ++st.fevals; return sys.deriv(tt, yy); };
    Vec<N> y1 = rk4_step(f, t, y, h);
    if (project) sys.project(y1);
    ++st.steps;
    return y1;
}

// Integrate from t0 until obs() returns false or t reaches tmax. obs(t, y) is called at t0 and
//  - for fixed-step methods: after every step (t = t0 + k*dt);
//  - for adaptive methods:   at the output grid t0 + k*dt using dense output, plus at tmax.
template <class Sys, class Obs>
IntegratorStats integrate(const Sys& sys, Method m, double t0, Vec<Sys::N> y, double tmax,
                          const IntegratorOptions& o, Obs&& obs) {
    constexpr std::size_t N = Sys::N;
    IntegratorStats st;
    double t = t0;
    if (!obs(t, y)) return st;
    if (!is_adaptive(m)) {
        long k = 0;
        while (t < tmax - 1e-12 && st.steps < o.max_steps) {
            y = fixed_step(sys, m, t, y, o.dt, st, o.project);
            ++k;
            t = t0 + k * o.dt;  // avoid accumulating round-off in t
            if (!obs(t, y)) break;
        }
        return st;
    }

    // ---- adaptive Dormand-Prince (plain or RKMK) ----
    if (is_lie(m) && !Sys::has_lie)
        throw std::invalid_argument("Lie-group integrators need a system with a rotational state (3D mode)");
    double h = (o.h0 > 0.0) ? o.h0 : o.dt;
    const double hmax = (o.hmax > 0.0) ? o.hmax : std::numeric_limits<double>::infinity();
    long kout = 1;
    double tout = t0 + o.dt, tlast_out = t0;
    Vec<N> k1{};
    bool have_k1 = false;
    while (t < tmax - 1e-12) {
        if (st.steps + st.rejected >= o.max_steps) break;
        h = std::min({h, hmax, tmax - t});
        double err = 0.0;
        bool stop = false;
        Vec<N> ynew;
        if (!is_lie(m)) {
            auto f = [&](double tt, const Vec<N>& yy) { ++st.fevals; return sys.deriv(tt, yy); };
            if (!have_k1) { k1 = f(t, y); have_k1 = true; }
            auto s = dopri_step(f, t, y, k1, h, o.rtol, o.atol);
            err = s.err;
            if (err <= 1.0 || h <= o.hmin) {
                while (tout <= t + h + 1e-12 && !stop) {
                    Vec<N> yi = s.dense((tout - t) / h);
                    if (o.project) sys.project(yi);
                    if (!obs(tout, yi)) stop = true;
                    tlast_out = tout;
                    tout = t0 + (++kout) * o.dt;
                }
                ynew = s.y1;
                if (o.project) sys.project(ynew);
                k1 = s.k7;  // FSAL; re-evaluate only if the projection actually moved the state
                if (ynew != s.y1) k1 = f(t + h, ynew);
            }
        } else {
            if constexpr (Sys::has_lie) {
                constexpr std::size_t NL = Sys::NL;
                const Vec<N> base = y;
                auto fl = [&](double tt, const Vec<NL>& z) { ++st.fevals; return sys.lie_deriv(tt, base, z); };
                const Vec<NL> z0 = sys.lie_zero(base);
                auto s = dopri_step(fl, t, z0, fl(t, z0), h, o.rtol, o.atol);
                err = s.err;
                if (err <= 1.0 || h <= o.hmin) {
                    while (tout <= t + h + 1e-12 && !stop) {
                        if (!obs(tout, sys.lie_to_global(base, s.dense((tout - t) / h)))) stop = true;
                        tlast_out = tout;
                        tout = t0 + (++kout) * o.dt;
                    }
                    ynew = sys.lie_to_global(base, s.y1);
                }
            }
        }
        if (err <= 1.0 || h <= o.hmin) {
            ++st.steps;
            if (stop) return st;
            t += h;
            y = ynew;
            if (t >= tmax - 1e-12 && t - tlast_out > 1e-9) obs(t, y);  // final point if off-grid
        } else {
            ++st.rejected;
        }
        // standard step-size controller (exponent 1/5, safety 0.9, growth limited to [0.2, 5])
        double fac = (err == 0.0) ? 5.0 : 0.9 * std::pow(err, -0.2);
        fac = std::min(5.0, std::max(0.2, fac));
        if (err > 1.0) fac = std::min(fac, 1.0);
        h = std::max(o.hmin, h * fac);
    }
    return st;
}

}  // namespace mk
