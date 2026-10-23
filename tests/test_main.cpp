// Self-contained test suite (no external framework):  make test
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "../src/integrators.hpp"
#include "../src/model.hpp"
#include "../src/sim2d.hpp"
#include "../src/sim3d.hpp"

using namespace mk;

static int g_fail = 0, g_pass = 0;
#define CHECK(cond, ...)                                                       \
    do {                                                                       \
        if (cond) { ++g_pass; }                                                \
        else { ++g_fail; std::printf("  FAIL %s:%d: %s  ", __FILE__, __LINE__, #cond); \
               std::printf(__VA_ARGS__); std::printf("\n"); }                  \
    } while (0)

// ---- helpers ------------------------------------------------------------------------------
template <std::size_t NN>
struct NoLie {
    static constexpr bool has_lie = false;
    static constexpr std::size_t NL = 1;
    void project(Vec<NN>&) const {}
    Vec<1> lie_zero(const Vec<NN>&) const { return {0.0}; }
    Vec<NN> lie_to_global(const Vec<NN>& b, const Vec<1>&) const { return b; }
    Vec<1> lie_deriv(double, const Vec<NN>&, const Vec<1>&) const { return {0.0}; }
};
// y' = y cos t, y(0) = 1  ->  y = exp(sin t)   (non-autonomous: checks the stage times)
struct ScalarODE : NoLie<1> {
    static constexpr std::size_t N = 1;
    Vec<1> deriv(double t, const Vec<1>& y) const { return {y[0] * std::cos(t)}; }
};
// harmonic oscillator x'' = -x
struct Oscillator : NoLie<2> {
    static constexpr std::size_t N = 2;
    Vec<2> deriv(double, const Vec<2>& y) const { return {y[1], -y[0]}; }
};

template <class Sys>
Vec<Sys::N> run_to(const Sys& sys, Method m, const Vec<Sys::N>& y0, double T, IntegratorOptions o,
                   IntegratorStats* st = nullptr) {
    Vec<Sys::N> last = y0;
    auto s = integrate(sys, m, 0.0, y0, T, o, [&](double, const Vec<Sys::N>& y) { last = y; return true; });
    if (st) *st = s;
    return last;
}
IntegratorOptions fixed(double dt) { IntegratorOptions o; o.dt = dt; return o; }
IntegratorOptions adaptive(double tol, double T) {
    IntegratorOptions o; o.rtol = o.atol = tol; o.dt = T; o.h0 = 1e-2; return o;
}
template <std::size_t N>
double maxabs(const Vec<N>& a, const Vec<N>& b) {
    double m = 0; for (std::size_t i = 0; i < N; ++i) m = std::max(m, std::fabs(a[i] - b[i])); return m;
}
// scaled error, matching the controller's norm
template <std::size_t N>
double scaled_err(const Vec<N>& a, const Vec<N>& ref) {
    double m = 0; for (std::size_t i = 0; i < N; ++i) m = std::max(m, std::fabs(a[i] - ref[i]) / (1.0 + std::fabs(ref[i])));
    return m;
}
double qangle(const Quat& a, const Quat& b) {  // rotation angle between two attitudes
    Quat d = qmul(qconj(a), b);
    return 2.0 * std::atan2(std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z), std::fabs(d.w));
}
Params3D torque_free(bool sphere) {
    Params3D p; p.aero = p.propulsion = p.gravity = false;
    if (sphere) { p.Ixx = p.Iyy = p.Izz = 1000.0; }
    else { p.Ixx = 948.0; p.Iyy = 1346.0; p.Izz = 1967.0; p.Ixz = 120.0; }
    return p;
}

// ---- tests ----------------------------------------------------------------------------------
void test_rk4_convergence() {
    std::printf("RK4 4th-order convergence\n");
    ScalarODE s; const double T = 5.0, exact = std::exp(std::sin(T));
    double prev = 0;
    for (double dt : {0.1, 0.05, 0.025, 0.0125}) {
        double e = std::fabs(run_to(s, Method::RK4, {1.0}, T, fixed(dt))[0] - exact);
        if (prev > 0) {
            double order = std::log2(prev / e);
            std::printf("  scalar  dt=%-7g err=%.3e observed order=%.3f\n", dt, e, order);
            CHECK(order > 3.8 && order < 4.2, "order %.3f", order);
        }
        prev = e;
    }
    // the actual 2D aircraft model (reference: rk4 with dt = 1/1280 s)
    Sys2D sys(1.5); const double T2 = 40.0;
    auto ref = run_to(sys, Method::RK4, Sys2D::initial(), T2, fixed(1.0 / 1280));
    prev = 0;
    for (double dt : {0.4, 0.2, 0.1, 0.05}) {
        double e = maxabs(run_to(sys, Method::RK4, Sys2D::initial(), T2, fixed(dt)), ref);
        if (prev > 0) {
            double order = std::log2(prev / e);
            std::printf("  2D model dt=%-6g err=%.3e observed order=%.3f\n", dt, e, order);
            CHECK(order > 3.7 && order < 4.4, "order %.3f", order);
        }
        prev = e;
    }
}

void test_2d_matches_original() {
    std::printf("2D refactor reproduces the original program (CL = 1.5, dt = 0.1)\n");
    Sys2D sys(1.5); IntegratorOptions o = fixed(0.1);
    Vec<4> y{}; double tf = 0;
    integrate(sys, Method::RK4, 0.0, Sys2D::initial(), 600.0, o, [&](double t, const Vec<4>& yy) {
        y = yy; tf = t; return yy[3] <= X_END; });
    // values printed by the original Matekane.cpp (commit b0d1f7f)
    std::printf("  V=%.6g gamma=%.6g h=%.6g x=%.6g t=%.4g\n", y[0], y[1], y[2], y[3], tf);
    CHECK(std::fabs(y[0] - 99.0293) < 1e-4 && std::fabs(y[1] - 0.0121954) < 1e-7 &&
          std::fabs(y[2] - 7727.87) < 1e-2 && std::fabs(y[3] - 10001.1) < 1e-1, "mismatch");
    CHECK(sweep_count(1.0, 2.0, 0.1) == 11 && sweep_count(1.0, 1.3, 0.1) == 4 && sweep_count(1, 2, 0) == 0,
          "sweep_count");
}

void test_dopri_coefficients() {
    std::printf("Dormand-Prince tableau: 5th-order fixed-step, dense output accuracy\n");
    ScalarODE s; double prev = 0, prevd = 0;
    for (int n : {40, 80, 160}) {
        double h = 5.0 / n, t = 0; Vec<1> y{1.0}; double dmax = 0;
        auto f = [&](double tt, const Vec<1>& yy) { return s.deriv(tt, yy); };
        for (int i = 0; i < n; ++i) {
            auto st = dopri_step(f, t, y, f(t, y), h, 1e-6, 1e-6);
            for (double th : {0.25, 0.5, 0.75}) {  // local dense error (start from exact solution)
                Vec<1> ye{std::exp(std::sin(t))};
                auto se = dopri_step(f, t, ye, f(t, ye), h, 1e-6, 1e-6);
                dmax = std::max(dmax, std::fabs(se.dense(th)[0] - std::exp(std::sin(t + th * h))));
            }
            y = st.y1; t += h;
        }
        double e = std::fabs(y[0] - std::exp(std::sin(5.0)));
        if (prev > 0) {
            double order = std::log2(prev / e), dorder = std::log2(prevd / dmax);
            std::printf("  h=%-8g global err=%.3e order=%.2f | dense local err=%.3e order=%.2f\n", h, e, order, dmax, dorder);
            CHECK(order > 4.6 && order < 5.6, "global order %.2f", order);
            CHECK(dorder > 4.3, "dense local order %.2f", dorder);  // local error O(h^5)
        }
        prev = e; prevd = dmax;
    }
}

void test_rk45_tolerance_and_efficiency() {
    std::printf("RK45: meets tolerance vs fine-step reference, fewer steps than RK4 for equal accuracy\n");
    Sys2D sys(1.5); const double T = 60.0;
    auto ref = run_to(sys, Method::RK4, Sys2D::initial(), T, fixed(1.0 / 2560));
    double prevErr = 1e9;
    for (double tol : {1e-6, 1e-8, 1e-10}) {
        IntegratorStats st;
        auto y = run_to(sys, Method::RK45, Sys2D::initial(), T, adaptive(tol, T), &st);
        double e = scaled_err(y, ref);
        std::printf("  tol=%-6g scaled err=%.3e  steps=%ld rejected=%ld fevals=%ld\n", tol, e, st.steps, st.rejected, st.fevals);
        CHECK(e < 50 * tol, "err %.3e > 50*tol", e);  // global error stays proportional to the local tolerance
        CHECK(e < prevErr, "error did not decrease with tolerance");
        prevErr = e;
    }
    // equal accuracy comparison at tol = 1e-8
    IntegratorStats s45;
    double e45 = scaled_err(run_to(sys, Method::RK45, Sys2D::initial(), T, adaptive(1e-8, T), &s45), ref);
    long n4 = 0, f4 = 0; double e4 = 1e9;
    for (int k = 0; k < 14 && e4 > e45; ++k) {
        double dt = 0.4 / std::pow(2.0, k);
        IntegratorStats s4;
        e4 = scaled_err(run_to(sys, Method::RK4, Sys2D::initial(), T, fixed(dt), &s4), ref);
        n4 = s4.steps; f4 = s4.fevals;
    }
    std::printf("  equal accuracy: rk45 err=%.2e with %ld steps/%ld fevals; rk4 needs %ld steps/%ld fevals (err=%.2e)\n",
                e45, s45.steps, s45.fevals, n4, f4, e4);
    CHECK(s45.steps < n4 && s45.fevals < f4, "rk45 not cheaper");
    // dense output on an oscillator: interpolated grid values meet tolerance
    Oscillator osc; double maxe = 0;
    IntegratorOptions o = adaptive(1e-9, 0.05); o.h0 = 0.1;
    integrate(osc, Method::RK45, 0.0, Vec<2>{1.0, 0.0}, 20.0, o, [&](double t, const Vec<2>& y) {
        maxe = std::max(maxe, std::fabs(y[0] - std::cos(t))); return true; });
    std::printf("  dense output (oscillator, tol 1e-9, output every 0.05 s): max err=%.2e\n", maxe);
    CHECK(maxe < 1e-7, "dense err %.2e", maxe);
}

void test_quaternion_normalization() {
    std::printf("Quaternion normalisation (rk4/rk45 renormalise each step; lie-* need none)\n");
    Params3D p; p.da = 3.0 * M_PI / 180; p.dr = -1.0 * M_PI / 180;  // banked, turning flight
    for (Method m : {Method::RK4, Method::RK45, Method::LieRK4, Method::LieRK45}) {
        Sys3D sys(p); double worst = 0;
        IntegratorOptions o = fixed(0.1); o.rtol = o.atol = 1e-9; o.h0 = 0.01;
        integrate(sys, m, 0.0, Sys3D::initial(), 60.0, o, [&](double, const Vec<13>& y) {
            worst = std::max(worst, std::fabs(qnorm(Sys3D::quat(y)) - 1.0)); return true; });
        std::printf("  %-8s max | |q|-1 | over 60 s rolling flight = %.2e\n", method_name(m), worst);
        CHECK(worst < 1e-14, "%s drift %.2e", method_name(m), worst);
    }
}

void test_pure_rotation() {
    std::printf("Pure rotation: constant body rate gives q(T) = q0 (x) exp(w T)\n");
    const V3 w{0.3, -0.2, 0.5};
    const Quat q0 = qfrom_euler(0.1, 0.2, 0.3);
    Sys3D sys(torque_free(true));  // spherical inertia: w stays constant
    const double T = 10.0;
    const Quat qex = qmul(q0, qexp({w[0] * T, w[1] * T, w[2] * T}));
    for (Method m : {Method::RK4, Method::RK45, Method::LieRK4, Method::LieRK45}) {
        IntegratorOptions o = fixed(0.01); o.rtol = o.atol = 1e-11; o.h0 = 0.01; o.dt = (is_adaptive(m) ? T : 0.01);
        auto y = run_to(sys, m, Sys3D::pack({0, 0, 0}, {0, 0, 0}, q0, w), T, o);
        double ang = qangle(Sys3D::quat(y), qex);
        std::printf("  %-8s attitude error after %.0f s = %.2e rad\n", method_name(m), T, ang);
        CHECK(ang < (is_lie(m) ? 1e-12 : 1e-9), "%s err %.2e", method_name(m), ang);
    }
    // yaw-only rotation: psi(T) = r T (wrapped)
    auto y = run_to(sys, Method::RK4, Sys3D::pack({0, 0, 0}, {0, 0, 0}, Quat{}, {0, 0, 0.25}), 4.0, fixed(0.05));
    V3 e = qto_euler(Sys3D::quat(y));
    std::printf("  yaw rate 0.25 rad/s for 4 s: psi=%.12f (expected 1), phi=%.1e theta=%.1e\n", e[2], e[0], e[1]);
    CHECK(std::fabs(e[2] - 1.0) < 1e-10 && std::fabs(e[0]) < 1e-14 && std::fabs(e[1]) < 1e-14, "yaw");
    // and the DCM convention: body x-axis rotated by +90 deg yaw points East in NED
    V3 xb = body_to_ned(qfrom_euler(0, 0, M_PI / 2), {1, 0, 0});
    CHECK(std::fabs(xb[1] - 1.0) < 1e-15 && std::fabs(xb[0]) < 1e-15, "dcm convention");
}

void test_lie_geometric() {
    std::printf("Geometric (RKMK) integrator: |q| = 1 without renormalisation, 4th order, invariants\n");
    Params3D p = torque_free(false); p.renormalize = false;
    Sys3D sys(p);
    const Vec<13> y0 = Sys3D::pack({0, 0, 0}, {0, 0, 0}, qfrom_euler(0.3, -0.4, 1.0), {0.05, 2.0, 0.05});  // tumbling (Dzhanibekov)
    const int nsteps = 100000; const double dt = 0.01;
    double worstLie = 0, worstRk4 = 0;
    integrate(sys, Method::LieRK4, 0.0, y0, nsteps * dt, fixed(dt), [&](double, const Vec<13>& y) {
        worstLie = std::max(worstLie, std::fabs(qnorm(Sys3D::quat(y)) - 1.0)); return true; });
    integrate(sys, Method::RK4, 0.0, y0, nsteps * dt, fixed(dt), [&](double, const Vec<13>& y) {
        worstRk4 = std::max(worstRk4, std::fabs(qnorm(Sys3D::quat(y)) - 1.0)); return true; });
    std::printf("  %d steps tumbling, NO renormalisation: lie-rk4 max| |q|-1 | = %.2e ; plain rk4 = %.2e\n",
                nsteps, worstLie, worstRk4);
    CHECK(worstLie < 1e-13, "lie drift %.2e", worstLie);
    CHECK(worstRk4 > worstLie, "expected plain rk4 to drift more");
    // lie-rk45 without renormalisation as well
    double worst45 = 0; IntegratorOptions o = adaptive(1e-10, 0.1); o.dt = 0.1;
    integrate(sys, Method::LieRK45, 0.0, y0, 500.0, o, [&](double, const Vec<13>& y) {
        worst45 = std::max(worst45, std::fabs(qnorm(Sys3D::quat(y)) - 1.0)); return true; });
    std::printf("  lie-rk45 500 s tumbling: max| |q|-1 | = %.2e\n", worst45);
    CHECK(worst45 < 1e-13, "lie45 drift %.2e", worst45);

    // convergence order of lie-rk4 against a fine reference (checks the dexp^-1 correction)
    const double T = 10.0;
    Params3D pr = p; pr.renormalize = true;
    auto ref = run_to(Sys3D(pr), Method::RK4, y0, T, fixed(1.0 / 2560));
    double prev = 0;
    for (double h : {0.08, 0.04, 0.02}) {
        auto y = run_to(sys, Method::LieRK4, y0, T, fixed(h));
        double e = qangle(Sys3D::quat(y), Sys3D::quat(ref)) + norm(Sys3D::omega(y) - Sys3D::omega(ref));
        if (prev > 0) {
            double order = std::log2(prev / e);
            std::printf("  lie-rk4 dt=%-5g err=%.3e order=%.2f\n", h, e, order);
            CHECK(order > 3.7 && order < 4.5, "order %.2f", order);
        }
        prev = e;
    }
    // lie-rk45 meets tolerance vs the reference
    auto y45 = run_to(sys, Method::LieRK45, y0, T, adaptive(1e-9, T));
    double e45 = qangle(Sys3D::quat(y45), Sys3D::quat(ref)) + norm(Sys3D::omega(y45) - Sys3D::omega(ref));
    std::printf("  lie-rk45 tol 1e-9: err=%.2e\n", e45);
    CHECK(e45 < 1e-7, "lie-rk45 err %.2e", e45);
    // torque-free invariants: inertial angular momentum vector and rotational energy
    auto H = [&](const Vec<13>& y) { V3 w = Sys3D::omega(y); V3 Iw{dot(sys.I[0], w), dot(sys.I[1], w), dot(sys.I[2], w)};
                                     return body_to_ned(Sys3D::quat(y), Iw); };
    auto E = [&](const Vec<13>& y) { V3 w = Sys3D::omega(y); return 0.5 * dot(w, V3{dot(sys.I[0], w), dot(sys.I[1], w), dot(sys.I[2], w)}); };
    auto y = run_to(sys, Method::LieRK4, y0, 100.0, fixed(0.01));
    double dH = norm(H(y) - H(y0)) / norm(H(y0)), dE = std::fabs(E(y) - E(y0)) / E(y0);
    std::printf("  lie-rk4 100 s: rel. change |H_ned| %.2e, energy %.2e\n", dH, dE);
    CHECK(dH < 1e-8 && dE < 1e-8, "invariants");
}

void test_ballistic() {
    std::printf("Free fall / ballistic (aero and thrust off)\n");
    for (bool spin : {false, true}) {
        Params3D p; p.aero = p.propulsion = false; p.Ixz = 120.0;
        Sys3D sys(p);
        const Quat q0 = qfrom_euler(0.2, 0.3, -0.5);
        const V3 vb0{150.0, -10.0, 20.0}, p0{0.0, 0.0, -H0};
        const V3 w0 = spin ? V3{0.4, 0.7, -0.3} : V3{0, 0, 0};
        const V3 vn0 = body_to_ned(q0, vb0);
        const double T = 20.0;
        for (Method m : {Method::RK4, Method::LieRK4, Method::RK45}) {
            IntegratorOptions o = fixed(0.01); if (is_adaptive(m)) o = adaptive(1e-11, T);
            auto y = run_to(sys, m, Sys3D::pack(p0, vb0, q0, w0), T, o);
            V3 pex{p0[0] + vn0[0] * T, p0[1] + vn0[1] * T, p0[2] + vn0[2] * T + 0.5 * G * T * T};
            V3 vex{vn0[0], vn0[1], vn0[2] + G * T};
            V3 vn = body_to_ned(Sys3D::quat(y), Sys3D::vel(y));
            double ep = norm(Sys3D::pos(y) - pex), ev = norm(vn - vex);
            std::printf("  %-8s %-9s pos err %.2e ft, vel err %.2e ft/s (drop %.1f ft)\n", method_name(m),
                        spin ? "tumbling" : "no spin", ep, ev, 0.5 * G * T * T);
            // exact for the non-spinning case (quadratic solution); with spin the body-frame
            // equations are nonlinear, so only truncation error (relative ~1e-9) remains
            CHECK(ep < (spin ? 1e-4 : 1e-8) && ev < (spin ? 1e-5 : 1e-9), "ballistic");
        }
    }
}

void test_3d_reduces_to_2d() {
    std::printf("3D constrained to the pitch plane reproduces 2D\n");
    const double CL = 1.5, T = 60.0, dt = 0.01;
    Sys2D s2(CL);
    std::vector<Vec<4>> traj2;
    integrate(s2, Method::RK4, 0.0, Sys2D::initial(), T, fixed(dt), [&](double, const Vec<4>& y) { traj2.push_back(y); return true; });
    const double gdot0 = s2.deriv(0, Sys2D::initial())[1];
    double prevErr = 1e9;
    for (double Cma : {-50.0, -500.0}) {
        Params3D p; p.CL_cmd = CL; p.Cma = Cma; p.CLq = 0.0;  // stiff pitch axis keeps alpha ~ 0
        Sys3D s3(p);
        double eh = 0, ex = 0, eV = 0, eg = 0, lat = 0; std::size_t k = 0;
        integrate(s3, Method::RK4, 0.0, Sys3D::initial(gdot0), T, fixed(dt), [&](double, const Vec<13>& y) {
            const auto& z = traj2[k++];
            V3 vn = body_to_ned(Sys3D::quat(y), Sys3D::vel(y));
            double V = norm(vn), gam = std::asin(-vn[2] / V);
            eh = std::max(eh, std::fabs(-y[2] - z[2])); ex = std::max(ex, std::fabs(y[0] - z[3]));
            eV = std::max(eV, std::fabs(V - z[0])); eg = std::max(eg, std::fabs(gam - z[1]));
            V3 e = qto_euler(Sys3D::quat(y));
            lat = std::max({lat, std::fabs(y[1]), std::fabs(y[4]), std::fabs(y[10]), std::fabs(y[12]), std::fabs(e[0]), std::fabs(e[2])});
            return true; });
        std::printf("  Cma=%-5g max |dh|=%.3f ft |dx|=%.3f ft |dV|=%.4f ft/s |dgamma|=%.2e rad ; lateral states max %.1e\n",
                    Cma, eh, ex, eV, eg, lat);
        CHECK(lat == 0.0, "lateral states excited");
        // the residual is the small angle of attack needed to pitch the aircraft (O(1/|Cma|))
        const double sc = 50.0 / std::fabs(Cma);
        CHECK(eh < 1.0 * sc && ex < 1.0 * sc && eV < 0.2 * sc && eg < 3e-3 * sc, "3D != 2D");
        CHECK(eh < 0.2 * prevErr, "error did not shrink ~1/|Cma| with a stiffer pitch axis");
        prevErr = eh;
    }
}

int main() {
    test_rk4_convergence();
    test_2d_matches_original();
    test_dopri_coefficients();
    test_rk45_tolerance_and_efficiency();
    test_quaternion_normalization();
    test_pure_rotation();
    test_lie_geometric();
    test_ballistic();
    test_3d_reduces_to_2d();
    std::printf("\n%d checks passed, %d failed\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
