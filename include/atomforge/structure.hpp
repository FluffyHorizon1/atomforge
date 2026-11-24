// atomforge — periodic atomic structure + neighbour list (header-only, C++17)
#pragma once
#include "math.hpp"
#include <map>
#include <string>
#include <vector>
#include <cmath>

namespace atomforge {

struct Structure {
    Mat3 cell;                    // rows = lattice vectors (Angstrom)
    std::vector<Vec3> pos;        // Cartesian positions (Angstrom)
    std::vector<int> species;     // index into the potential's element list
    std::string label;            // e.g. "L1_2", "fcc", "mc-ordered"

    std::size_t size() const { return pos.size(); }
    double volume() const { return std::abs(det(cell)); }

    Vec3 to_frac(const Vec3& r) const { return transpose(inverse(cell)) * r; }
    Vec3 to_cart(const Vec3& f) const { return transpose(cell) * f; }

    void add_frac(const Vec3& f, int sp) { pos.push_back(to_cart(f)); species.push_back(sp); }

    void wrap() {
        const Mat3 inv_t = transpose(inverse(cell));
        const Mat3 cell_t = transpose(cell);
        for (auto& r : pos) {
            Vec3 f = inv_t * r;
            for (int k = 0; k < 3; ++k) f[k] -= std::floor(f[k]);
            r = cell_t * f;
        }
    }

    // Apply a deformation F (r -> F r) to cell and positions.
    void deform(const Mat3& F) {
        cell = transpose(F * transpose(cell));
        for (auto& r : pos) r = F * r;
    }

    Structure supercell(int n1, int n2, int n3) const {
        Structure s;
        s.label = label;
        s.cell = Mat3::rows(cell.row(0) * n1, cell.row(1) * n2, cell.row(2) * n3);
        for (int a = 0; a < n1; ++a) for (int b = 0; b < n2; ++b) for (int c = 0; c < n3; ++c) {
            const Vec3 shift = cell.row(0) * a + cell.row(1) * b + cell.row(2) * c;
            for (std::size_t i = 0; i < size(); ++i) { s.pos.push_back(pos[i] + shift); s.species.push_back(species[i]); }
        }
        return s;
    }

    std::map<int, int> composition() const {
        std::map<int, int> c; for (int s : species) ++c[s]; return c;
    }
    double fraction(int sp) const {
        int n = 0; for (int s : species) n += (s == sp); return size() ? double(n) / size() : 0.0;
    }
};

// Perpendicular widths of the cell (distance between opposite faces).
inline Vec3 face_widths(const Mat3& cell) {
    const Vec3 a = cell.row(0), b = cell.row(1), c = cell.row(2);
    const double V = std::abs(dot(a, cross(b, c)));
    return {V / norm(cross(b, c)), V / norm(cross(c, a)), V / norm(cross(a, b))};
}

// Half neighbour list (each pair once) over all periodic images within rc.
// Brute force O(N^2 * images): fine for the cell sizes used in structure search (<~1000 atoms).
struct NeighborList {
    struct Pair { int i, j; Vec3 shift; };   // r_ij = pos[j] + shift - pos[i]
    std::vector<Pair> pairs;
    double cutoff{0}, skin{0};
    std::vector<Vec3> ref_pos;
    Mat3 ref_cell;

    void build(const Structure& s, double rc, double skin_ = 0.3) {
        cutoff = rc; skin = skin_;
        pairs.clear();
        const double rl = rc + skin, rl2 = rl * rl;
        const Vec3 w = face_widths(s.cell);
        const int n0 = int(std::ceil(rl / w.x)), n1 = int(std::ceil(rl / w.y)), n2 = int(std::ceil(rl / w.z));
        const int N = int(s.size());
        for (int a = -n0; a <= n0; ++a) for (int b = -n1; b <= n1; ++b) for (int c = -n2; c <= n2; ++c) {
            const Vec3 shift = s.cell.row(0) * a + s.cell.row(1) * b + s.cell.row(2) * c;
            const bool zero = (a == 0 && b == 0 && c == 0);
            for (int i = 0; i < N; ++i) for (int j = (zero ? i + 1 : i); j < N; ++j) {
                // For i==j with nonzero shift, count only one of (shift, -shift).
                if (i == j && !(a > 0 || (a == 0 && (b > 0 || (b == 0 && c > 0))))) continue;
                const Vec3 d = s.pos[j] + shift - s.pos[i];
                if (dot(d, d) < rl2) pairs.push_back({i, j, shift});
            }
        }
        ref_pos = s.pos; ref_cell = s.cell;
    }

    // True if atoms moved / cell changed enough that the list may be stale.
    bool needs_rebuild(const Structure& s) const {
        if (ref_pos.size() != s.size()) return true;
        for (int k = 0; k < 9; ++k) if (std::abs(ref_cell.m[k] - s.cell.m[k]) > 1e-9) return true;
        double mx = 0;
        for (std::size_t i = 0; i < s.size(); ++i) { const Vec3 d = s.pos[i] - ref_pos[i]; mx = std::max(mx, dot(d, d)); }
        return std::sqrt(mx) > 0.5 * skin;
    }
};

} // namespace atomforge
