#ifndef CONJUNCTION_ASSESSMENT_H
#define CONJUNCTION_ASSESSMENT_H

/**
 * Conjunction Assessment Engine
 *
 * Pipeline: TLE pairs → SGP4 propagation → find TCA → B-plane projection
 *           → Alfano maximum collision probability → CDM output
 *
 * Collision probability method: Alfano (AAS 03-548)
 *   "Relating Position Uncertainty to Maximum Conjunction Probability"
 *   PDF: https://celestrak.org/SOCRATES/AIAA-03-548.pdf
 *
 * Validation target: CelesTrak SOCRATES Plus
 *   https://celestrak.org/SOCRATES/
 *   Method: SGP4 propagator, 5km threshold, 7-day lookahead
 *   Default RTN covariance: 100m R, 300m T, 100m N
 *
 * Covariance source:
 *   - Default: SOCRATES-style fixed RTN (100m, 300m, 100m)
 *   - OD from OEM: via Tudat's OrbitDeterminationManager + fitOrbitToEphemeris
 *     Source: DigitalArsenal/tudat-wasm
 *     Files: src/tudatpy_wasm/estimation/estimation_analysis/
 *   - Propagated: via Tudat's propagateCovarianceRsw()
 */

#include "sgp4_propagator.h"
#include <vector>
#include <string>

namespace conjunction {

/// Default hard-body radii (meters) — SOCRATES uses combined ~10m for LEO objects
constexpr double DEFAULT_RADIUS_M = 5.0;

/// Default RTN covariance (meters) — matches SOCRATES
constexpr double DEFAULT_COV_R_M = 100.0;  // Radial
constexpr double DEFAULT_COV_T_M = 300.0;  // In-track (along-track)
constexpr double DEFAULT_COV_N_M = 100.0;  // Cross-track (normal)

/// Screening threshold (km)
constexpr double DEFAULT_THRESHOLD_KM = 5.0;

/// Conjunction event
struct ConjunctionEvent {
    // Object identifiers
    TLE obj1;
    TLE obj2;

    // Time of closest approach
    double tca_jd = 0.0;
    std::string tca_iso;

    // At TCA
    double min_range_km = 0.0;
    double rel_speed_kms = 0.0;

    // States at TCA (TEME)
    StateVector state1;
    StateVector state2;

    // Relative position/velocity at TCA (RTN frame)
    double rel_pos_r = 0.0, rel_pos_t = 0.0, rel_pos_n = 0.0;
    double rel_vel_r = 0.0, rel_vel_t = 0.0, rel_vel_n = 0.0;

    // Collision probability (Alfano maximum probability method)
    double max_probability = 0.0;
    double dilution_threshold_km = 0.0;
    std::string probability_method = "ALFANO-MAXPROB";

    // Covariance (RTN, meters²)
    double cov_r1 = 0.0, cov_t1 = 0.0, cov_n1 = 0.0;  // Object 1
    double cov_r2 = 0.0, cov_t2 = 0.0, cov_n2 = 0.0;  // Object 2

    // Days since epoch for each object
    double dse1 = 0.0, dse2 = 0.0;
};

/// Find Time of Closest Approach between two TLEs
/// Searches within [start_jd, start_jd + duration_days]
/// Uses bisection refinement after coarse step search
double find_tca(const TLE& tle1, const TLE& tle2,
                double start_jd, double duration_days = 7.0,
                double coarse_step_sec = 60.0, double fine_tol_sec = 0.001);

/// Compute full conjunction assessment for a TLE pair
ConjunctionEvent assess_conjunction(
    const TLE& tle1, const TLE& tle2,
    double start_jd, double duration_days = 7.0,
    double radius1_m = DEFAULT_RADIUS_M,
    double radius2_m = DEFAULT_RADIUS_M);

/// Compute full conjunction assessment near a known TCA
/// Searches ±window_hours around tca_hint_jd
ConjunctionEvent assess_conjunction_near(
    const TLE& tle1, const TLE& tle2,
    double tca_hint_jd, double window_hours = 2.0,
    double radius1_m = DEFAULT_RADIUS_M,
    double radius2_m = DEFAULT_RADIUS_M);

/// Compute Alfano maximum collision probability
/// d = miss distance (km), Rc = combined hard-body radius (km)
/// Returns {max_probability, dilution_threshold_km}
struct ProbResult {
    double max_probability;
    double dilution_threshold_km;
    double sigma_star_km;
};
ProbResult alfano_max_probability(double miss_distance_km, double combined_radius_km);

/// Compute collision probability with covariance
/// Full B-plane projection method
/// pos1, vel1, pos2, vel2: states in same frame (km, km/s)
/// cov1, cov2: 3x3 position covariance matrices (km²)
/// combined_radius: hard-body radius (km)
double collision_probability(
    const StateVector& state1, const StateVector& state2,
    const double cov1[9], const double cov2[9],
    double combined_radius_km);

/// Screen a set of TLEs for conjunctions
std::vector<ConjunctionEvent> screen_conjunctions(
    const std::vector<TLE>& primary_tles,
    const std::vector<TLE>& secondary_tles,
    double start_jd,
    double duration_days = 7.0,
    double threshold_km = DEFAULT_THRESHOLD_KM);

/// Convert inertial state to RTN frame relative to a reference state
void inertial_to_rtn(const StateVector& ref, const StateVector& target,
                     double& r, double& t, double& n,
                     double& vr, double& vt, double& vn);

} // namespace conjunction

#endif // CONJUNCTION_ASSESSMENT_H
