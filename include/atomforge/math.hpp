// atomforge — small linear algebra for 3D atomistics (header-only, C++17)
#pragma once
#include <array>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace atomforge {

struct Vec3 {
    double x{0}, y{0}, z{0};
    constexpr Vec3() = default;
    constexpr Vec3(double a, double b, double c) : x(a), y(b), z(c) {}
    double& operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }
    double operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }
    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
};
inline Vec3 operator+(Vec3 a, const Vec3& b) { return a += b; }
inline Vec3 operator-(Vec3 a, const Vec3& b) { return a -= b; }
inline Vec3 operator-(const Vec3& a) { return {-a.x, -a.y, -a.z}; }
inline Vec3 operator*(Vec3 a, double s) { return a *= s; }
inline Vec3 operator*(double s, Vec3 a) { return a *= s; }
inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

// Row-major 3x3. When used as a cell, row i is lattice vector i.
struct Mat3 {
    std::array<double, 9> m{0, 0, 0, 0, 0, 0, 0, 0, 0};
    double& operator()(int i, int j) { return m[3 * i + j]; }
    double operator()(int i, int j) const { return m[3 * i + j]; }
    static Mat3 identity() { Mat3 r; r(0, 0) = r(1, 1) = r(2, 2) = 1; return r; }
    static Mat3 diag(double a, double b, double c) { Mat3 r; r(0, 0) = a; r(1, 1) = b; r(2, 2) = c; return r; }
    static Mat3 rows(const Vec3& a, const Vec3& b, const Vec3& c) {
        Mat3 r; for (int j = 0; j < 3; ++j) { r(0, j) = a[j]; r(1, j) = b[j]; r(2, j) = c[j]; } return r;
    }
    Vec3 row(int i) const { return {(*this)(i, 0), (*this)(i, 1), (*this)(i, 2)}; }
};
inline Mat3 operator*(const Mat3& a, const Mat3& b) {
    Mat3 r;
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
        double s = 0; for (int k = 0; k < 3; ++k) s += a(i, k) * b(k, j); r(i, j) = s;
    }
    return r;
}
inline Mat3 operator+(Mat3 a, const Mat3& b) { for (int i = 0; i < 9; ++i) a.m[i] += b.m[i]; return a; }
inline Mat3 operator-(Mat3 a, const Mat3& b) { for (int i = 0; i < 9; ++i) a.m[i] -= b.m[i]; return a; }
inline Mat3 operator*(Mat3 a, double s) { for (auto& v : a.m) v *= s; return a; }
inline Vec3 operator*(const Mat3& a, const Vec3& v) {
    return {a(0, 0) * v.x + a(0, 1) * v.y + a(0, 2) * v.z,
            a(1, 0) * v.x + a(1, 1) * v.y + a(1, 2) * v.z,
            a(2, 0) * v.x + a(2, 1) * v.y + a(2, 2) * v.z};
}
inline Mat3 transpose(const Mat3& a) { Mat3 r; for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) r(i, j) = a(j, i); return r; }
inline double det(const Mat3& a) {
    return a(0, 0) * (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)) - a(0, 1) * (a(1, 0) * a(2, 2) - a(1, 2) * a(2, 0)) +
           a(0, 2) * (a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0));
}
inline Mat3 inverse(const Mat3& a) {
    const double d = det(a);
    if (std::abs(d) < 1e-14) throw std::runtime_error("singular matrix");
    Mat3 r;
    r(0, 0) = (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)) / d;
    r(0, 1) = (a(0, 2) * a(2, 1) - a(0, 1) * a(2, 2)) / d;
    r(0, 2) = (a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1)) / d;
    r(1, 0) = (a(1, 2) * a(2, 0) - a(1, 0) * a(2, 2)) / d;
    r(1, 1) = (a(0, 0) * a(2, 2) - a(0, 2) * a(2, 0)) / d;
    r(1, 2) = (a(0, 2) * a(1, 0) - a(0, 0) * a(1, 2)) / d;
    r(2, 0) = (a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0)) / d;
    r(2, 1) = (a(0, 1) * a(2, 0) - a(0, 0) * a(2, 1)) / d;
    r(2, 2) = (a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0)) / d;
    return r;
}

// Symmetric eigenvalues by cyclic Jacobi (small N, used for 6x6 elastic tensors).
template <std::size_t N>
std::array<double, N> sym_eigenvalues(std::array<std::array<double, N>, N> a) {
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0;
        for (std::size_t p = 0; p < N; ++p) for (std::size_t q = p + 1; q < N; ++q) off += a[p][q] * a[p][q];
        if (off < 1e-22) break;
        for (std::size_t p = 0; p < N; ++p) for (std::size_t q = p + 1; q < N; ++q) {
            if (std::abs(a[p][q]) < 1e-300) continue;
            const double theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
            const double t = (theta >= 0 ? 1.0 : -1.0) / (std::abs(theta) + std::sqrt(theta * theta + 1));
            const double c = 1 / std::sqrt(t * t + 1), s = t * c;
            for (std::size_t k = 0; k < N; ++k) {
                const double akp = a[k][p], akq = a[k][q];
                a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq;
            }
            for (std::size_t k = 0; k < N; ++k) {
                const double apk = a[p][k], aqk = a[q][k];
                a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk;
            }
        }
    }
    std::array<double, N> ev{};
    for (std::size_t i = 0; i < N; ++i) ev[i] = a[i][i];
    std::sort(ev.begin(), ev.end());
    return ev;
}

namespace units {
inline constexpr double eV_per_A3_to_GPa = 160.21766208;
inline constexpr double kB = 8.617333262e-5;            // eV/K
inline constexpr double time_fs = 10.180505710774743;   // 1 sqrt(amu*A^2/eV) in fs
} // namespace units

} // namespace atomforge
