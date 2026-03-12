/**
 * Conjunction Assessment Engine (v2) — Propagator-agnostic implementation
 *
 * TCA finding uses the same golden-section refinement as v1, but
 * through the EphemerisSource interface instead of TLE+SGP4 directly.
 *
 * B-plane projection follows Montenbruck & Gill Section 6.5.
 */

#include "conjunction/conjunction_engine.h"
#include "conjunction/conjunction_assessment.h"  // for inertial_to_rtn

#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace conjunction {

// ── Constructor ──

ConjunctionEngine::ConjunctionEngine()
    : pc_method_(std::make_unique<Foster2D>()),
      default_cov_(Covariance3x3::socrates_default()) {}

void ConjunctionEngine::set_pc_method(std::unique_ptr<PcMethod> method) {
    pc_method_ = std::move(method);
}

void ConjunctionEngine::set_pc_method(const std::string& name) {
    pc_method_ = create_pc_method(name);
}

void ConjunctionEngine::set_combined_radius_m(double r1, double r2) {
    radius1_m_ = r1;
    radius2_m_ = r2;
}

void ConjunctionEngine::set_default_covariance(const Covariance3x3& cov) {
    default_cov_ = cov;
}

std::string ConjunctionEngine::pc_method_name() const {
    return pc_method_ ? pc_method_->name() : "NONE";
}

// ── Helpers ──

static double vec_norm(const double v[3]) {
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static void vec_normalize(double v[3]) {
    double n = vec_norm(v);
    if (n > 1e-15) { v[0] /= n; v[1] /= n; v[2] /= n; }
}

static double vec_dot(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static void vec_cross(const double a[3], const double b[3], double c[3]) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

// Matrix-vector multiply: out = M × v (3×3 × 3)
static void mat_vec(const double M[9], const double v[3], double out[3]) {
    out[0] = M[0]*v[0] + M[1]*v[1] + M[2]*v[2];
    out[1] = M[3]*v[0] + M[4]*v[1] + M[5]*v[2];
    out[2] = M[6]*v[0] + M[7]*v[1] + M[8]*v[2];
}

// R^T × C × R for 3×3 rotation and covariance → 2×2 projection
// R is 3×2 (two columns = encounter-plane basis vectors)

static double distance_at(const EphemerisSource& o1, const EphemerisSource& o2, double jd) {
    auto s1 = o1.state_at(jd);
    auto s2 = o2.state_at(jd);
    double dx = s1.x - s2.x, dy = s1.y - s2.y, dz = s1.z - s2.z;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

// ── TCA Finding ──

double ConjunctionEngine::find_tca(
    const EphemerisSource& obj1, const EphemerisSource& obj2,
    double start_jd, double duration_days, double coarse_step_sec) const
{
    double end_jd = start_jd + duration_days;
    double step = coarse_step_sec / 86400.0;

    // Phase 1: Find all local minima
    struct Minimum { double jd; double dist; };
    std::vector<Minimum> minima;

    double prev_d = 1e18, curr_d = 1e18, prev_jd = start_jd;
    bool prev_decreasing = false;

    for (double jd = start_jd; jd <= end_jd; jd += step) {
        try {
            double d = distance_at(obj1, obj2, jd);
            bool decreasing = (d < curr_d);
            if (prev_decreasing && !decreasing && curr_d < 50.0) {
                minima.push_back({prev_jd, curr_d});
            }
            prev_decreasing = decreasing;
            prev_d = curr_d;
            curr_d = d;
            prev_jd = jd;
        } catch (...) {}
    }

    if (minima.empty()) return start_jd;

    // Phase 2: Sort and refine top candidates
    std::sort(minima.begin(), minima.end(),
              [](const Minimum& a, const Minimum& b) { return a.dist < b.dist; });

    int n_refine = std::min(static_cast<int>(minima.size()), 30);
    double best_jd = minima[0].jd;
    double best_dist = 1e18;

    for (int i = 0; i < n_refine; i++) {
        double center = minima[i].jd;

        // Sub-second rescan ±10s at 0.05s
        double subscan_step = 0.05 / 86400.0;
        double subscan_window = 10.0 / 86400.0;
        double sub_best_jd = center, sub_best_d = minima[i].dist;

        for (double jd = center - subscan_window; jd <= center + subscan_window; jd += subscan_step) {
            try {
                double d = distance_at(obj1, obj2, jd);
                if (d < sub_best_d) { sub_best_d = d; sub_best_jd = jd; }
            } catch (...) {}
        }

        // Golden section refinement ±1s
        double a = sub_best_jd - 1.0/86400.0;
        double b = sub_best_jd + 1.0/86400.0;
        double tol = 0.001 / 86400.0;
        const double phi = (std::sqrt(5.0) - 1.0) / 2.0;

        while ((b - a) > tol) {
            double c = b - phi * (b - a);
            double d = a + phi * (b - a);
            try {
                double fc = distance_at(obj1, obj2, c);
                double fd = distance_at(obj1, obj2, d);
                if (fc < fd) b = d; else a = c;
            } catch (...) { break; }
        }

        double refined_jd = (a + b) / 2.0;
        try {
            double d = distance_at(obj1, obj2, refined_jd);
            if (d < best_dist) { best_dist = d; best_jd = refined_jd; }
        } catch (...) {}
    }

    return best_jd;
}

// ── B-plane Construction ──

BPlaneGeometry ConjunctionEngine::build_bplane(
    const StateVector& s1, const StateVector& s2,
    const Covariance3x3& c1, const Covariance3x3& c2,
    double combined_radius_km) const
{
    BPlaneGeometry bp;
    bp.combined_radius = combined_radius_km;

    // Relative position and velocity
    double dr[3] = {s1.x - s2.x, s1.y - s2.y, s1.z - s2.z};
    double dv[3] = {s1.vx - s2.vx, s1.vy - s2.vy, s1.vz - s2.vz};
    double v_rel = vec_norm(dv);

    if (v_rel < 1e-10) {
        bp.xi = vec_norm(dr);
        bp.sigma_xx = 1e-6;
        bp.sigma_zz = 1e-6;
        return bp;
    }

    // Encounter-plane basis (⊥ to relative velocity)
    double zhat[3] = {dv[0]/v_rel, dv[1]/v_rel, dv[2]/v_rel};

    // Find xhat perpendicular to zhat
    double temp[3] = {1, 0, 0};
    if (std::abs(zhat[0]) > 0.9) { temp[0] = 0; temp[1] = 1; }

    double d_tz = vec_dot(temp, zhat);
    double xhat[3] = {temp[0] - d_tz*zhat[0], temp[1] - d_tz*zhat[1], temp[2] - d_tz*zhat[2]};
    vec_normalize(xhat);

    double yhat[3];
    vec_cross(zhat, xhat, yhat);

    // Project miss vector onto encounter plane
    bp.xi = vec_dot(dr, xhat);
    bp.zeta = vec_dot(dr, yhat);

    // Combined covariance
    double C[9];
    for (int i = 0; i < 9; i++) C[i] = c1.data[i] + c2.data[i];

    // Project covariance onto encounter plane (2×2)
    double Cx[3], Cy[3];
    mat_vec(C, xhat, Cx);
    mat_vec(C, yhat, Cy);

    bp.sigma_xx = vec_dot(xhat, Cx);
    bp.sigma_xz = vec_dot(xhat, Cy);
    bp.sigma_zz = vec_dot(yhat, Cy);

    return bp;
}

// ── Full Assessment ──

ConjunctionEvent2 ConjunctionEngine::assess(
    const EphemerisSource& obj1, const EphemerisSource& obj2,
    double start_jd, double duration_days,
    const Covariance3x3* cov1, const Covariance3x3* cov2) const
{
    ConjunctionEvent2 event;

    // Metadata
    event.obj1_name = obj1.object_name();
    event.obj2_name = obj2.object_name();
    event.obj1_id = obj1.object_id();
    event.obj2_id = obj2.object_id();
    event.obj1_norad = obj1.norad_id();
    event.obj2_norad = obj2.norad_id();

    // Find TCA
    event.tca_jd = find_tca(obj1, obj2, start_jd, duration_days);
    event.tca_iso = jd_to_iso(event.tca_jd);

    // Get states at TCA
    event.state1 = obj1.state_at(event.tca_jd);
    event.state2 = obj2.state_at(event.tca_jd);

    // Geometry
    double dx = event.state1.x - event.state2.x;
    double dy = event.state1.y - event.state2.y;
    double dz = event.state1.z - event.state2.z;
    event.miss_distance_km = std::sqrt(dx*dx + dy*dy + dz*dz);

    double dvx = event.state1.vx - event.state2.vx;
    double dvy = event.state1.vy - event.state2.vy;
    double dvz = event.state1.vz - event.state2.vz;
    event.relative_speed_kms = std::sqrt(dvx*dvx + dvy*dvy + dvz*dvz);

    // RTN
    inertial_to_rtn(event.state1, event.state2,
                    event.rel_r, event.rel_t, event.rel_n,
                    event.rel_vr, event.rel_vt, event.rel_vn);

    // DSE
    event.dse1 = event.tca_jd - obj1.epoch_jd();
    event.dse2 = event.tca_jd - obj2.epoch_jd();

    // Covariance
    event.cov1 = cov1 ? *cov1 : default_cov_;
    event.cov2 = cov2 ? *cov2 : default_cov_;
    event.combined_radius_km = (radius1_m_ + radius2_m_) / 1000.0;

    // Build B-plane and compute Pc
    event.bplane = build_bplane(event.state1, event.state2,
                                 event.cov1, event.cov2,
                                 event.combined_radius_km);
    event.pc = pc_method_->compute(event.bplane);

    return event;
}

ConjunctionEvent2 ConjunctionEngine::assess_near(
    const EphemerisSource& obj1, const EphemerisSource& obj2,
    double tca_hint_jd, double window_hours,
    const Covariance3x3* cov1, const Covariance3x3* cov2) const
{
    double window_days = window_hours / 24.0;
    return assess(obj1, obj2, tca_hint_jd - window_days, 2.0 * window_days,
                  cov1, cov2);
}

ConjunctionEvent2 ConjunctionEngine::compute_pc(
    const StateVector& state1, const StateVector& state2,
    const Covariance3x3& cov1, const Covariance3x3& cov2,
    double combined_radius_km) const
{
    ConjunctionEvent2 event;
    event.state1 = state1;
    event.state2 = state2;
    event.tca_jd = state1.epoch_jd;
    event.tca_iso = jd_to_iso(event.tca_jd);
    event.cov1 = cov1;
    event.cov2 = cov2;
    event.combined_radius_km = combined_radius_km;

    double dx = state1.x - state2.x;
    double dy = state1.y - state2.y;
    double dz = state1.z - state2.z;
    event.miss_distance_km = std::sqrt(dx*dx + dy*dy + dz*dz);

    double dvx = state1.vx - state2.vx;
    double dvy = state1.vy - state2.vy;
    double dvz = state1.vz - state2.vz;
    event.relative_speed_kms = std::sqrt(dvx*dvx + dvy*dvy + dvz*dvz);

    inertial_to_rtn(state1, state2,
                    event.rel_r, event.rel_t, event.rel_n,
                    event.rel_vr, event.rel_vt, event.rel_vn);

    event.bplane = build_bplane(state1, state2, cov1, cov2, combined_radius_km);
    event.pc = pc_method_->compute(event.bplane);

    return event;
}

std::vector<ConjunctionEvent2> ConjunctionEngine::screen(
    const std::vector<std::shared_ptr<EphemerisSource>>& primaries,
    const std::vector<std::shared_ptr<EphemerisSource>>& secondaries,
    double start_jd, double duration_days, double threshold_km) const
{
    std::vector<ConjunctionEvent2> events;

    for (const auto& primary : primaries) {
        for (const auto& secondary : secondaries) {
            if (primary->norad_id() != 0 && primary->norad_id() == secondary->norad_id())
                continue;

            try {
                auto event = assess(*primary, *secondary, start_jd, duration_days);
                if (event.miss_distance_km <= threshold_km) {
                    events.push_back(event);
                }
            } catch (...) {}
        }
    }

    std::sort(events.begin(), events.end(),
              [](const ConjunctionEvent2& a, const ConjunctionEvent2& b) {
                  return a.pc.probability > b.pc.probability;
              });

    return events;
}

} // namespace conjunction
