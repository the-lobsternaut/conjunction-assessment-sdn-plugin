#ifndef CONJUNCTION_PC_METHOD_H
#define CONJUNCTION_PC_METHOD_H

/**
 * PcMethod — Pluggable Collision Probability Computation
 *
 * All methods receive the same B-plane geometry and produce Pc.
 * The conjunction assessor is method-agnostic.
 *
 * Implementations:
 *   AlfanoMaxPc       — Maximum probability bound (AAS 03-548)
 *   Foster2D          — Foster-Estes short-encounter (1992), 2D Gaussian
 *   Patera2001        — Patera single-integral method (2001)
 *   Chan2008          — Chan series expansion (2008)
 *   Alfriend2D        — Alfriend-Akella 2D numerical integration
 *
 * References:
 *   Alfano: "Relating Position Uncertainty to Maximum Conjunction Probability"
 *           AAS 03-548, celestrak.org/SOCRATES/AIAA-03-548.pdf
 *   Foster: "A Parametric Analysis of Orbital Debris Collision Probability
 *            and Maneuver Rate for Space Vehicles" NASA/JSC-25898, 1992
 *   Patera: "General Method for Calculating Satellite Collision Probability"
 *            JGCD Vol 24 No 4, 2001
 *   Chan:   "Spacecraft Collision Probability" El Segundo: Aerospace Press, 2008
 *   Alfriend: "Probability of Collision Error Analysis" Space Debris Vol 1, 1999
 */

#include <string>
#include <cmath>

namespace conjunction {

// ── B-plane Geometry (input to all Pc methods) ──

struct BPlaneGeometry {
    // Miss distance components in encounter plane (km)
    double xi = 0;        // miss along encounter-plane X
    double zeta = 0;      // miss along encounter-plane Z

    // Combined covariance in encounter plane (km²)
    // 2×2 symmetric: [sxx, sxz; sxz, szz]
    double sigma_xx = 0;
    double sigma_xz = 0;
    double sigma_zz = 0;

    // Combined hard-body radius (km)
    double combined_radius = 0.01;  // 10m default

    // Miss distance magnitude (km)
    double miss_distance() const { return std::sqrt(xi * xi + zeta * zeta); }

    // ── Mahalanobis Distance (2D, encounter plane) ──
    // Md = sqrt(Δr^T × C⁻¹ × Δr)  in the B-plane
    // Tells you how many "sigma" the miss is from the covariance center
    double mahalanobis_distance() const {
        double det = sigma_xx * sigma_zz - sigma_xz * sigma_xz;
        if (det < 1e-30) return 0;
        // C⁻¹ = (1/det) × [σ_zz, -σ_xz; -σ_xz, σ_xx]
        double inv_xx = sigma_zz / det;
        double inv_xz = -sigma_xz / det;
        double inv_zz = sigma_xx / det;
        double md2 = xi * xi * inv_xx + 2.0 * xi * zeta * inv_xz + zeta * zeta * inv_zz;
        return (md2 > 0) ? std::sqrt(md2) : 0;
    }

    // Mahalanobis distance along principal axes
    double mahalanobis_principal() const {
        double s1, s2;
        double l1, l2;
        eigenvalues(l1, l2);
        s1 = (l1 > 0) ? l1 : 1e-30;
        s2 = (l2 > 0) ? l2 : 1e-30;
        double u, v;
        double theta = rotation_angle();
        double ct = std::cos(theta), st = std::sin(theta);
        u = ct * xi + st * zeta;
        v = -st * xi + ct * zeta;
        return std::sqrt(u * u / s1 + v * v / s2);
    }

    // Covariance eigenvalues
    void eigenvalues(double& lambda1, double& lambda2) const {
        double trace = sigma_xx + sigma_zz;
        double det = sigma_xx * sigma_zz - sigma_xz * sigma_xz;
        double disc = std::sqrt(std::max(0.0, trace * trace / 4.0 - det));
        lambda1 = trace / 2.0 + disc;
        lambda2 = trace / 2.0 - disc;
    }

    // Eigenvalue standard deviations
    void eigen_sigmas(double& s1, double& s2) const {
        double l1, l2;
        eigenvalues(l1, l2);
        s1 = (l1 > 0) ? std::sqrt(l1) : 0;
        s2 = (l2 > 0) ? std::sqrt(l2) : 0;
    }

    // Rotation angle to diagonalize covariance
    double rotation_angle() const {
        if (std::abs(sigma_xz) < 1e-30 && std::abs(sigma_xx - sigma_zz) < 1e-30)
            return 0;
        return 0.5 * std::atan2(2.0 * sigma_xz, sigma_xx - sigma_zz);
    }

    // Rotated miss vector (along principal axes)
    void rotated_miss(double& u, double& v) const {
        double theta = rotation_angle();
        double ct = std::cos(theta), st = std::sin(theta);
        u = ct * xi + st * zeta;
        v = -st * xi + ct * zeta;
    }
};

// ── Pc Result ──

struct PcResult {
    double probability = 0;
    std::string method;
    bool converged = true;
    int iterations = 0;          // for iterative methods
    double max_probability = 0;  // Alfano upper bound (always computed)
    double mahalanobis_2d = 0;   // Mahalanobis distance in encounter plane
};

// ── Abstract Interface ──

class PcMethod {
public:
    virtual ~PcMethod() = default;

    /// Compute collision probability from B-plane geometry
    virtual PcResult compute(const BPlaneGeometry& bplane) const = 0;

    /// Method name
    virtual std::string name() const = 0;
};

// ── Alfano Maximum Probability (AAS 03-548) ──

class AlfanoMaxPc : public PcMethod {
public:
    PcResult compute(const BPlaneGeometry& bplane) const override;
    std::string name() const override { return "ALFANO-MAXPROB"; }
};

// ── Foster-Estes 2D (1992) ──
// Short-encounter approximation: 2D Gaussian integration over hard-body disk.
// Industry standard for CDM Pc values. Accurate when encounter duration
// is short relative to covariance evolution.
//
// P_c = (R²/2σ₁σ₂) × exp(-½(u²/σ₁² + v²/σ₂²))
// where u,v are miss components along covariance eigenvectors

class Foster2D : public PcMethod {
public:
    PcResult compute(const BPlaneGeometry& bplane) const override;
    std::string name() const override { return "FOSTER-2D"; }
};

// ── Patera 2001 (single integral) ──
// Reduces the 2D probability integral to a 1D integral along the
// miss-distance azimuth. More accurate than Foster for highly
// elliptical covariance.
//
// Ref: Patera, "General Method for Calculating Satellite Collision Probability"
// JGCD 24(4), 2001

class Patera2001 : public PcMethod {
public:
    explicit Patera2001(int quadrature_points = 64) : n_quad_(quadrature_points) {}
    PcResult compute(const BPlaneGeometry& bplane) const override;
    std::string name() const override { return "PATERA-2001"; }
private:
    int n_quad_;
};

// ── Chan 2008 (series expansion) ──
// Exact series solution using modified Bessel functions.
// Converges fast for moderate aspect ratios. The standard
// for NASA CARA operations.
//
// Ref: Chan, "Spacecraft Collision Probability" Aerospace Press, 2008

class Chan2008 : public PcMethod {
public:
    explicit Chan2008(int max_terms = 100, double tol = 1e-16) : max_terms_(max_terms), tol_(tol) {}
    PcResult compute(const BPlaneGeometry& bplane) const override;
    std::string name() const override { return "CHAN-2008"; }
private:
    int max_terms_;
    double tol_;
};

// ── Alfriend-Akella 2D Numerical Integration ──
// Direct numerical integration of the 2D Gaussian PDF over the
// hard-body disk. Most general method — works for any covariance
// shape. Uses Gauss-Legendre quadrature.
//
// Ref: Alfriend & Akella, "Probability of Collision Error Analysis"
// Space Debris 1(1), 1999

class Alfriend2D : public PcMethod {
public:
    explicit Alfriend2D(int radial_points = 32, int angular_points = 64)
        : n_r_(radial_points), n_theta_(angular_points) {}
    PcResult compute(const BPlaneGeometry& bplane) const override;
    std::string name() const override { return "ALFRIEND-2D"; }
private:
    int n_r_, n_theta_;
};

// ── Factory ──

/// Create Pc method by name
/// Names: "alfano", "foster", "patera", "chan", "alfriend"
std::unique_ptr<PcMethod> create_pc_method(const std::string& name);

} // namespace conjunction

#endif // CONJUNCTION_PC_METHOD_H
