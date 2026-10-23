/* Matekane Airstrip takeoff simulation
** created December 8, 2011 by metafalsevac
**
** 2D mode : original vertical-plane point-mass model (wind axes, flight-path angle).
** 3D mode : 6-DOF rigid-body model with unit-quaternion attitude (src/sim3d.hpp).
** Integrators: rk4 (default), rk45 (adaptive Dormand-Prince), lie-rk4 / lie-rk45 (3D only,
** Runge-Kutta-Munthe-Kaas: quaternion stays exactly on the unit sphere).
**
** Run `matekane --help` for options. With no --cl* options the program prompts for the
** lift-coefficient sweep interactively, like the original.
*/
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include "src/integrators.hpp"
#include "src/model.hpp"
#include "src/sim2d.hpp"
#include "src/sim3d.hpp"

using namespace mk;

namespace {

constexpr double DEG = M_PI / 180.0;

struct Config {
    std::string mode = "2d";
    Method method = Method::RK4;
    IntegratorOptions io;
    double tmax = 600.0;  // safety limit on simulated time, s
    bool have_cl = false;
    double cl_min = 0, cl_max = 0, cl_step = 0.1;
    std::string out;
    Params3D p3;
    bool stats = false;
};

void usage() {
    std::cout <<
R"(Usage: matekane [options]
  --mode 2d|3d             simulation model (default 2d)
  --integrator NAME        rk4 (default) | rk45 | lie-rk4 | lie-rk45   (lie-* need --mode 3d)
  --dt S                   rk4 step and output interval, s (default 0.1)
  --rtol X --atol X        rk45 tolerances (default 1e-8, 1e-8)
  --cl X                   simulate a single lift coefficient
  --cl-min X --cl-max X --cl-step X   lift coefficient sweep (prompted for if omitted)
  --tmax S                 simulated-time limit per run, s (default 600)
  --out FILE               trajectory output (default output.txt for 2d, output3d.csv for 3d)
  3D only:
  --aileron DEG --elevator DEG --rudder DEG   constant control deflections
  --no-renorm              do not renormalise the quaternion after rk4/rk45 steps
  --stats                  print integrator step / function-evaluation counts
)";
}

double num(int& i, int argc, char** argv) {
    if (i + 1 >= argc) throw std::invalid_argument(std::string("missing value for ") + argv[i]);
    char* end = nullptr;
    double v = std::strtod(argv[++i], &end);
    if (!end || *end) throw std::invalid_argument(std::string("bad number: ") + argv[i]);
    return v;
}

Config parse(int argc, char** argv) {
    Config c;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-h" || a == "--help") { usage(); std::exit(0); }
        else if (a == "--mode" && i + 1 < argc) c.mode = argv[++i];
        else if (a == "--integrator" && i + 1 < argc) c.method = parse_method(argv[++i]);
        else if (a == "--dt") c.io.dt = num(i, argc, argv);
        else if (a == "--rtol") c.io.rtol = num(i, argc, argv);
        else if (a == "--atol") c.io.atol = num(i, argc, argv);
        else if (a == "--tmax") c.tmax = num(i, argc, argv);
        else if (a == "--cl") { c.cl_min = c.cl_max = num(i, argc, argv); c.cl_step = 1.0; c.have_cl = true; }
        else if (a == "--cl-min") { c.cl_min = num(i, argc, argv); c.have_cl = true; }
        else if (a == "--cl-max") { c.cl_max = num(i, argc, argv); c.have_cl = true; }
        else if (a == "--cl-step") c.cl_step = num(i, argc, argv);
        else if (a == "--out" && i + 1 < argc) c.out = argv[++i];
        else if (a == "--aileron") c.p3.da = num(i, argc, argv) * DEG;
        else if (a == "--elevator") c.p3.de = num(i, argc, argv) * DEG;
        else if (a == "--rudder") c.p3.dr = num(i, argc, argv) * DEG;
        else if (a == "--no-renorm") c.p3.renormalize = false;
        else if (a == "--stats") c.stats = true;
        else throw std::invalid_argument("unknown option " + a + " (see --help)");
    }
    if (c.mode != "2d" && c.mode != "3d") throw std::invalid_argument("--mode must be 2d or 3d");
    if (c.mode == "2d" && is_lie(c.method))
        throw std::invalid_argument("lie-rk4/lie-rk45 act on the quaternion attitude and require --mode 3d");
    if (!(c.io.dt > 0) || !(c.tmax > 0) || !(c.io.rtol > 0) || !(c.io.atol > 0))
        throw std::invalid_argument("dt, tmax, rtol and atol must be positive");
    if (c.out.empty()) c.out = (c.mode == "2d") ? "output.txt" : "output3d.csv";
    return c;
}

enum class Status { Crash, Survive, NoResult };

template <class Sys>
struct Outcome {
    Status status = Status::NoResult;
    double t = 0.0;
    double x_touch = 0.0;
    Vec<Sys::N> y{};
    IntegratorStats st;
};

// Integrate one case until the aircraft hits the terrain, passes X_END, loses airspeed or tmax.
template <class Sys, class Writer>
Outcome<Sys> run_case(const Sys& sys, const Vec<Sys::N>& y0, const Config& c, Writer&& write) {
    Outcome<Sys> out;
    double xPrev = Sys::dist(y0), cPrev = Sys::alt(y0) - ground(Sys::dist(y0));
    bool first = true;
    auto obs = [&](double t, const Vec<Sys::N>& y) {
        write(t, y);
        out.t = t;
        out.y = y;
        bool finite = true;
        for (double v : y) finite = finite && std::isfinite(v);
        if (!finite || Sys::airspeed(y) <= 0.0) { out.status = Status::NoResult; return false; }
        const double x = Sys::dist(y), clr = Sys::alt(y) - ground(x);
        if (!first && clr <= 0.0) {  // crashed: interpolate the touchdown point between samples
            out.status = Status::Crash;
            out.x_touch = xPrev + (x - xPrev) * (cPrev / (cPrev - clr));
            return false;
        }
        if (x > X_END) { out.status = Status::Survive; return false; }
        first = false;
        xPrev = x;
        cPrev = clr;
        return true;
    };
    out.st = integrate(sys, c.method, 0.0, y0, c.tmax, c.io, obs);
    return out;
}

void print_stats(const Config& c, const IntegratorStats& st) {
    if (!c.stats) return;
    std::cout << "[" << method_name(c.method) << ": " << st.steps << " steps, " << st.rejected << " rejected, "
              << st.fevals << " f-evals]" << std::endl;
}

const char* status_text(Status s) {
    return s == Status::Crash ? " CRASH!!!" : (s == Status::Survive ? " SURVIVE!!!" : " NO RESULT (airspeed lost or time limit reached)");
}

void run2d(const Config& c, double CL, std::ofstream& f) {
    Sys2D sys(CL);
    auto w = [&](double, const Vec<4>& y) {
        f << y[3] << "     " << y[2] << "     " << y[0] << "      " << ground(y[3]) << "     " << y[1] << "     "
          << veas(y[0], y[2]) << "     " << CL << "\n";
    };
    auto o = run_case(sys, Sys2D::initial(), c, w);
    const auto& y = o.y;
    std::cout << "CL = " << CL << status_text(o.status) << std::endl;
    std::cout << "Vxw = " << y[0] << " gamma = " << y[1] << " Alt = " << y[2] << " x = " << y[3]
              << " Elev = " << ground(y[3]) << "\nEquiv Airspeed = " << veas(y[0], y[2]);
    if (o.status == Status::Crash) std::cout << "\nTouchdown (interpolated) x = " << o.x_touch;
    std::cout << "  t = " << o.t << " s" << std::endl;
    print_stats(c, o.st);
    std::cout << std::endl;
}

void run3d(const Config& c, double CL, std::ofstream& f) {
    Params3D p = c.p3;
    p.CL_cmd = CL;
    Sys3D sys(p);
    auto w = [&](double t, const Vec<13>& y) {
        const Quat q = Sys3D::quat(y);
        const V3 e = qto_euler(q);
        const AeroState a = sys.aero_state(y);
        f << t << ',' << y[0] << ',' << y[1] << ',' << y[2] << ',' << -y[2] << ',' << ground(y[0]);
        for (int i = 3; i < 13; ++i) f << ',' << y[i];
        f << ',' << e[0] / DEG << ',' << e[1] / DEG << ',' << e[2] / DEG << ',' << a.alpha / DEG << ','
          << a.beta / DEG << ',' << a.V << ',' << a.CL << ',' << qnorm(q) - 1.0 << ',' << CL << "\n";
    };
    auto o = run_case(sys, Sys3D::initial(), c, w);
    const auto& y = o.y;
    const V3 e = qto_euler(Sys3D::quat(y));
    const V3 vn = body_to_ned(Sys3D::quat(y), Sys3D::vel(y));
    const double V = norm(vn), gam = std::asin(-vn[2] / V);
    std::cout << "CL = " << CL << status_text(o.status) << std::endl;
    std::cout << "V = " << V << " gamma = " << gam << " Alt = " << -y[2] << " x = " << y[0] << " y = " << y[1]
              << " Elev = " << ground(y[0]) << "\nroll/pitch/yaw (deg) = " << e[0] / DEG << " / " << e[1] / DEG
              << " / " << e[2] / DEG << "  alpha (deg) = " << sys.aero_state(y).alpha / DEG
              << "  |q|-1 = " << qnorm(Sys3D::quat(y)) - 1.0;
    if (o.status == Status::Crash) std::cout << "\nTouchdown (interpolated) x = " << o.x_touch;
    std::cout << "  t = " << o.t << " s" << std::endl;
    print_stats(c, o.st);
    std::cout << std::endl;
}

}  // namespace

int main(int argc, char** argv) {
    Config c;
    try {
        c = parse(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << std::endl;
        return 2;
    }

    std::cout << "Takeoff From Matekane Airstrip (" << c.mode << " mode, " << method_name(c.method) << ")\n\n";
    if (!c.have_cl) {
        std::cout << "This program will simulate takeoff from the Matekane Airstrip \nat various lift coefficients.\n\n"
                     "Input minimum lift coefficient: ";
        std::cin >> c.cl_min;
        std::cout << "\nInput maximum lift coefficient: ";
        std::cin >> c.cl_max;
        std::cout << "\nInput lift coefficient step size: ";
        std::cin >> c.cl_step;
        std::cout << "\n" << std::endl;
        if (!std::cin) { std::cerr << "error: invalid input" << std::endl; return 2; }
    }
    if (!(c.cl_step > 0.0) || c.cl_max < c.cl_min) {
        std::cerr << "error: need CLmin <= CLmax and step > 0" << std::endl;
        return 2;
    }

    std::ofstream f(c.out);
    if (c.mode == "2d") {
        f << "# x(ft)  h(ft)  V(ft/s)  ground(ft)  gamma(rad)  Veas(ft/s)  CL\n";
    } else {
        f << "t,x_N,y_E,z_D,h,ground,u,v,w,qw,qx,qy,qz,p,q,r,phi_deg,theta_deg,psi_deg,alpha_deg,beta_deg,V,CL,"
             "qnorm_err,CL_cmd\n";
    }
    f << std::setprecision(10);

    const int n = sweep_count(c.cl_min, c.cl_max, c.cl_step);
    for (int i = 0; i < n; ++i) {
        const double CL = c.cl_min + i * c.cl_step;
        if (CL > CL_MAX + 1e-12) {
            std::cout << "CL = " << CL << " exceeds CLMax = " << CL_MAX << " (stalled wing) - skipped\n" << std::endl;
            continue;
        }
        if (c.mode == "2d") run2d(c, CL, f);
        else run3d(c, CL, f);
    }
    std::cout << "Trajectory written to " << c.out << std::endl;
    return 0;
}
