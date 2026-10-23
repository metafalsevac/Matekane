// Small fixed-size vector and quaternion helpers (header-only, no dependencies).
#pragma once
#include <array>
#include <cmath>
#include <cstddef>

namespace mk {

template <std::size_t N>
using Vec = std::array<double, N>;
using V3 = Vec<3>;

template <std::size_t N>
inline Vec<N> operator+(const Vec<N>& a, const Vec<N>& b) {
    Vec<N> r;
    for (std::size_t i = 0; i < N; ++i) r[i] = a[i] + b[i];
    return r;
}
template <std::size_t N>
inline Vec<N> operator-(const Vec<N>& a, const Vec<N>& b) {
    Vec<N> r;
    for (std::size_t i = 0; i < N; ++i) r[i] = a[i] - b[i];
    return r;
}
template <std::size_t N>
inline Vec<N> operator*(double s, const Vec<N>& a) {
    Vec<N> r;
    for (std::size_t i = 0; i < N; ++i) r[i] = s * a[i];
    return r;
}
template <std::size_t N>
inline Vec<N> operator*(const Vec<N>& a, double s) { return s * a; }

template <std::size_t N>
inline double dot(const Vec<N>& a, const Vec<N>& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < N; ++i) s += a[i] * b[i];
    return s;
}
template <std::size_t N>
inline double norm(const Vec<N>& a) { return std::sqrt(dot(a, a)); }

inline V3 cross(const V3& a, const V3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

// ---------------------------------------------------------------------------
// Quaternions. Convention (documented in README):
//   * Hamilton algebra (i*j = k), SCALAR-FIRST storage q = [w, x, y, z].
//   * q is the attitude of the BODY frame relative to the inertial NED frame and
//     rotates body-frame vectors into NED:  v_ned = q (x) [0, v_body] (x) q*  =  R(q) v_body.
//   * Kinematics with body angular rate w_b:  q_dot = 0.5 * q (x) [0, w_b].
// ---------------------------------------------------------------------------
struct Quat {
    double w = 1.0, x = 0.0, y = 0.0, z = 0.0;
};

inline Quat qmul(const Quat& a, const Quat& b) {
    return {a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
            a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
            a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
            a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w};
}
inline Quat qconj(const Quat& q) { return {q.w, -q.x, -q.y, -q.z}; }
inline double qnorm(const Quat& q) { return std::sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z); }
inline Quat qnormalized(const Quat& q) {
    double n = qnorm(q);
    return {q.w / n, q.x / n, q.y / n, q.z / n};
}

// Rotation matrix body -> NED (rows are NED axes expressed in body components).
inline std::array<V3, 3> qdcm(const Quat& q) {
    const double w = q.w, x = q.x, y = q.y, z = q.z;
    return {{{1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)},
             {2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)},
             {2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)}}};
}
inline V3 body_to_ned(const Quat& q, const V3& vb) {
    auto R = qdcm(q);
    return {dot(R[0], vb), dot(R[1], vb), dot(R[2], vb)};
}
inline V3 ned_to_body(const Quat& q, const V3& vn) {
    auto R = qdcm(q);  // R^T * vn
    return {R[0][0] * vn[0] + R[1][0] * vn[1] + R[2][0] * vn[2],
            R[0][1] * vn[0] + R[1][1] * vn[1] + R[2][1] * vn[2],
            R[0][2] * vn[0] + R[1][2] * vn[1] + R[2][2] * vn[2]};
}

// Quaternion exponential of a rotation vector theta (rad): [cos(|t|/2), sin(|t|/2) t/|t|].
// Always returns a unit quaternion (to rounding).
inline Quat qexp(const V3& th) {
    double a = norm(th);
    double h = 0.5 * a;
    double s = (a < 1e-8) ? 0.5 - a * a / 48.0 : std::sin(h) / a;  // sin(a/2)/a
    return {std::cos(h), s * th[0], s * th[1], s * th[2]};
}

// 3-2-1 (yaw psi, pitch theta, roll phi) Euler angles <-> quaternion.
inline Quat qfrom_euler(double phi, double theta, double psi) {
    double cr = std::cos(phi / 2), sr = std::sin(phi / 2);
    double cp = std::cos(theta / 2), sp = std::sin(theta / 2);
    double cy = std::cos(psi / 2), sy = std::sin(psi / 2);
    return {cr * cp * cy + sr * sp * sy, sr * cp * cy - cr * sp * sy, cr * sp * cy + sr * cp * sy,
            cr * cp * sy - sr * sp * cy};
}
inline V3 qto_euler(const Quat& q) {  // returns {phi, theta, psi}
    double sinp = 2 * (q.w * q.y - q.z * q.x);
    sinp = sinp > 1 ? 1 : (sinp < -1 ? -1 : sinp);
    return {std::atan2(2 * (q.w * q.x + q.y * q.z), 1 - 2 * (q.x * q.x + q.y * q.y)), std::asin(sinp),
            std::atan2(2 * (q.w * q.z + q.x * q.y), 1 - 2 * (q.y * q.y + q.z * q.z))};
}

// Inverse of the right Jacobian of SO(3): if R(t) = R0 * Exp(theta(t)) and R^T R_dot = [w]x
// (body rates), then theta_dot = Jr^{-1}(theta) w
//   = w + 1/2 theta x w + c(|theta|) theta x (theta x w),  c = (1 - (a/2) cot(a/2)) / a^2.
inline V3 dexpinv_right(const V3& th, const V3& w) {
    double a2 = dot(th, th);
    double c;
    if (a2 < 1e-6) {
        c = 1.0 / 12.0 + a2 / 720.0 + a2 * a2 / 30240.0;
    } else {
        double a = std::sqrt(a2);
        c = (1.0 - 0.5 * a / std::tan(0.5 * a)) / a2;
    }
    V3 txw = cross(th, w);
    V3 txtxw = cross(th, txw);
    return {w[0] + 0.5 * txw[0] + c * txtxw[0], w[1] + 0.5 * txw[1] + c * txtxw[1],
            w[2] + 0.5 * txw[2] + c * txtxw[2]};
}

}  // namespace mk
