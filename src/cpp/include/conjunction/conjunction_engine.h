#ifndef CONJUNCTION_ENGINE_H
#define CONJUNCTION_ENGINE_H

/**
 * Conjunction Assessment Engine (v2) — Propagator & Pc-method agnostic
 *
 * Decoupled architecture:
 *   EphemerisSource → provides states at any time
 *   PcMethod        → computes collision probability
 *   ConjunctionEngine → orchestrates TCA finding + Pc computation
 *
 * The engine never knows or cares what propagator produced the states
 * or what method computed the probability.
 */

#include "conjunction/ephemeris_source.h"
#include "conjunction/pc_method.h"
#include <memory>
#include <vector>

namespace conjunction {

// ── Covariance ──

/// 3×3 position covariance matrix (km²), row-major
struct Covariance3x3 {
    double data[9] = {0};

    /// Construct from RTN diagonal (σ_R, σ_T, σ_N in km)
    static Covariance3x3 from_rtn_diagonal(double sr, double st, double sn) {
        Covariance3x3 c;
        c.data[0] = sr * sr; c.data[4] = st * st; c.data[8] = sn * sn;
        return c;
    }

    /// Construct from full 3×3 matrix
    static Covariance3x3 from_matrix(const double m[9]) {
        Covariance3x3 c;
        for (int i = 0; i < 9; i++) c.data[i] = m[i];
        return c;
    }

    /// Default SOCRATES covariance (100m R, 300m T, 100m N)
    static Covariance3x3 socrates_default() {
        return from_rtn_diagonal(0.1, 0.3, 0.1);
    }
};

// ── Conjunction Event (v2) ──

struct ConjunctionEvent2 {
    // Object identifiers
    std::string obj1_name, obj2_name;
    std::string obj1_id, obj2_id;
    int obj1_norad = 0, obj2_norad = 0;

    // Time of closest approach
    double tca_jd = 0;
    std::string tca_iso;

    // States at TCA (inertial frame)
    StateVector state1, state2;

    // Geometry
    double miss_distance_km = 0;
    double relative_speed_kms = 0;

    // RTN relative position/velocity
    double rel_r = 0, rel_t = 0, rel_n = 0;
    double rel_vr = 0, rel_vt = 0, rel_vn = 0;

    // B-plane geometry
    BPlaneGeometry bplane;

    // Collision probability (from chosen method)
    PcResult pc;

    // Covariance used
    Covariance3x3 cov1, cov2;
    double combined_radius_km = 0.01;

    // Days since epoch
    double dse1 = 0, dse2 = 0;
};

// ── Engine ──

class ConjunctionEngine {
public:
    ConjunctionEngine();

    /// Set the Pc method (default: Foster-2D)
    void set_pc_method(std::unique_ptr<PcMethod> method);
    void set_pc_method(const std::string& name);

    /// Set default hard-body radii (meters)
    void set_combined_radius_m(double radius1, double radius2);

    /// Set default covariance
    void set_default_covariance(const Covariance3x3& cov);

    /// Find TCA between two ephemeris sources
    /// Searches [start_jd, start_jd + duration_days]
    double find_tca(const EphemerisSource& obj1, const EphemerisSource& obj2,
                     double start_jd, double duration_days = 7.0,
                     double coarse_step_sec = 5.0) const;

    /// Full conjunction assessment
    ConjunctionEvent2 assess(
        const EphemerisSource& obj1, const EphemerisSource& obj2,
        double start_jd, double duration_days = 7.0,
        const Covariance3x3* cov1 = nullptr,
        const Covariance3x3* cov2 = nullptr) const;

    /// Assess near a known TCA
    ConjunctionEvent2 assess_near(
        const EphemerisSource& obj1, const EphemerisSource& obj2,
        double tca_hint_jd, double window_hours = 2.0,
        const Covariance3x3* cov1 = nullptr,
        const Covariance3x3* cov2 = nullptr) const;

    /// Compute Pc from pre-computed states + covariance (no propagation)
    ConjunctionEvent2 compute_pc(
        const StateVector& state1, const StateVector& state2,
        const Covariance3x3& cov1, const Covariance3x3& cov2,
        double combined_radius_km = 0.01) const;

    /// Screen multiple objects
    std::vector<ConjunctionEvent2> screen(
        const std::vector<std::shared_ptr<EphemerisSource>>& primaries,
        const std::vector<std::shared_ptr<EphemerisSource>>& secondaries,
        double start_jd, double duration_days = 7.0,
        double threshold_km = 5.0) const;

    /// Get current Pc method name
    std::string pc_method_name() const;

private:
    std::unique_ptr<PcMethod> pc_method_;
    double radius1_m_ = 5.0;
    double radius2_m_ = 5.0;
    Covariance3x3 default_cov_;

    /// Build B-plane geometry from states + covariance
    BPlaneGeometry build_bplane(
        const StateVector& s1, const StateVector& s2,
        const Covariance3x3& c1, const Covariance3x3& c2,
        double combined_radius_km) const;
};

} // namespace conjunction

#endif // CONJUNCTION_ENGINE_H
