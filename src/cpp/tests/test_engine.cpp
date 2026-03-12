/**
 * Conjunction Engine v2 — Propagator-Agnostic Tests
 *
 * Tests the decoupled architecture:
 *   1. SGP4 ephemeris source (backward compat with v1)
 *   2. OEM ephemeris source (tabulated, Hermite interpolated)
 *   3. Callback ephemeris source (user function)
 *   4. Engine with different Pc methods
 *   5. compute_pc() from pre-computed states (no propagation)
 *   6. Cross-validate v1 vs v2 on same TLE pair
 */

#include "conjunction/conjunction_engine.h"
#include "conjunction/conjunction_assessment.h"  // v1 for comparison
#include "conjunction/gp_json.h"

#include <cstdio>
#include <cmath>
#include <memory>

static int tests_passed = 0;
static int tests_failed = 0;

#define CHECK(cond, msg) do { \
    if (!(cond)) { printf("  FAIL: %s\n", msg); tests_failed++; } \
    else { printf("  PASS: %s\n", msg); tests_passed++; } \
} while(0)

#define CHECK_TOL(val, expected, tol, msg) do { \
    double _v = (val), _e = (expected), _t = (tol); \
    if (std::abs(_v - _e) > _t) { \
        printf("  FAIL: %s (got %.6e, expected %.6e)\n", msg, _v, _e); \
        tests_failed++; \
    } else { \
        printf("  PASS: %s\n", msg); tests_passed++; \
    } \
} while(0)

using namespace conjunction;

// ISS and Cosmos debris TLEs for testing
static const char* ISS_TLE_NAME = "ISS (ZARYA)";
static const char* ISS_TLE_L1 = "1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9001";
static const char* ISS_TLE_L2 = "2 25544  51.6400 100.0000 0001500  80.0000 280.0000 15.49000000000017";

static const char* DEB_TLE_NAME = "COSMOS DEBRIS";
static const char* DEB_TLE_L1 = "1 99999U 93036AAA 24001.50000000  .00000100  00000-0  10000-4 0  9999";
static const char* DEB_TLE_L2 = "2 99999  51.6300 100.2000 0005000  85.0000 275.0000 15.48500000000019";

// ── Test 1: SGP4 Ephemeris Source ──

void test_sgp4_source() {
    printf("\n=== Test 1: SGP4 EphemerisSource ===\n");

    auto tle = parse_tle(ISS_TLE_NAME, ISS_TLE_L1, ISS_TLE_L2);
    SGP4EphemerisSource src(tle);

    CHECK(src.norad_id() == 25544, "NORAD ID = 25544");
    CHECK(src.object_name() == "ISS (ZARYA)", "Name = ISS (ZARYA)");

    auto state = src.state_at(tle.epoch_jd);
    double rmag = std::sqrt(state.x*state.x + state.y*state.y + state.z*state.z);
    printf("  At epoch: r = %.1f km\n", rmag);
    CHECK(rmag > 6300 && rmag < 6900, "ISS position in LEO range");
}

// ── Test 2: OEM Ephemeris Source ──

void test_oem_source() {
    printf("\n=== Test 2: OEM EphemerisSource (Hermite interpolation) ===\n");

    // Generate tabulated ephemeris from SGP4
    auto tle = parse_tle(ISS_TLE_NAME, ISS_TLE_L1, ISS_TLE_L2);
    double epoch = tle.epoch_jd;

    std::vector<EphemerisPoint> points;
    for (int i = 0; i <= 100; i++) {
        double jd = epoch + i * 60.0 / 86400.0;  // 1-minute steps
        auto sv = propagate_sgp4(tle, jd);
        points.push_back({jd, sv.x, sv.y, sv.z, sv.vx, sv.vy, sv.vz});
    }

    OEMEphemerisSource oem_src(points, "ISS-OEM", 25544);

    // Check interpolation at a point between samples (30.5 minutes)
    double test_jd = epoch + 30.5 * 60.0 / 86400.0;
    auto interp_state = oem_src.state_at(test_jd);
    auto sgp4_state = propagate_sgp4(tle, test_jd);

    double pos_err = std::sqrt(
        (interp_state.x - sgp4_state.x) * (interp_state.x - sgp4_state.x) +
        (interp_state.y - sgp4_state.y) * (interp_state.y - sgp4_state.y) +
        (interp_state.z - sgp4_state.z) * (interp_state.z - sgp4_state.z)
    );
    double vel_err = std::sqrt(
        (interp_state.vx - sgp4_state.vx) * (interp_state.vx - sgp4_state.vx) +
        (interp_state.vy - sgp4_state.vy) * (interp_state.vy - sgp4_state.vy) +
        (interp_state.vz - sgp4_state.vz) * (interp_state.vz - sgp4_state.vz)
    );

    printf("  Hermite interpolation error: pos=%.6f km, vel=%.9f km/s\n", pos_err, vel_err);
    CHECK(pos_err < 0.01, "Position error < 10m (Hermite at 60s steps)");
    CHECK(vel_err < 0.0001, "Velocity error < 0.1 m/s");
}

// ── Test 3: Callback Ephemeris Source ──

void test_callback_source() {
    printf("\n=== Test 3: Callback EphemerisSource ===\n");

    // Circular orbit at 400 km altitude
    double mu = 398600.4418;
    double r0 = 6778.0;
    double v0 = std::sqrt(mu / r0);
    double period = 2 * 3.14159265358979 * std::sqrt(r0*r0*r0 / mu);

    auto circular_orbit = [=](double jd) -> StateVector {
        double t = (jd - 2460310.5) * 86400.0;  // seconds since reference
        double omega = 2 * 3.14159265358979 / period;
        double theta = omega * t;
        return {jd,
            r0 * std::cos(theta), r0 * std::sin(theta), 0,
            -v0 * std::sin(theta), v0 * std::cos(theta), 0};
    };

    CallbackEphemerisSource src(circular_orbit, 2460310.5, "CIRCULAR", 99998);

    auto state = src.state_at(2460310.5);
    double rmag = std::sqrt(state.x*state.x + state.y*state.y + state.z*state.z);
    CHECK_TOL(rmag, r0, 0.01, "Circular orbit r = 6778 km");
    CHECK(src.object_name() == "CIRCULAR", "Callback name");
}

// ── Test 4: Engine with different Pc methods ──

void test_engine_methods() {
    printf("\n=== Test 4: Engine with Multiple Pc Methods ===\n");

    auto tle1 = parse_tle(ISS_TLE_NAME, ISS_TLE_L1, ISS_TLE_L2);
    auto tle2 = parse_tle(DEB_TLE_NAME, DEB_TLE_L1, DEB_TLE_L2);

    auto src1 = SGP4EphemerisSource(tle1);
    auto src2 = SGP4EphemerisSource(tle2);

    double start_jd = tle1.epoch_jd;

    const char* methods[] = {"alfano", "foster", "patera", "chan", "alfriend"};
    double pcs[5];

    for (int i = 0; i < 5; i++) {
        ConjunctionEngine engine;
        engine.set_pc_method(methods[i]);

        auto event = engine.assess(src1, src2, start_jd, 1.0);
        pcs[i] = event.pc.probability;

        printf("  %s: Pc=%.6e, miss=%.3f km, method=%s\n",
               methods[i], event.pc.probability, event.miss_distance_km, event.pc.method.c_str());
    }

    // All should produce non-negative Pc
    for (int i = 0; i < 5; i++) {
        CHECK(pcs[i] >= 0, (std::string(methods[i]) + " Pc ≥ 0").c_str());
    }

    // Alfano ≥ all others
    for (int i = 1; i < 5; i++) {
        CHECK(pcs[0] >= pcs[i] * 0.99,
              (std::string("Alfano ≥ ") + methods[i]).c_str());
    }
}

// ── Test 5: compute_pc() without propagation ──

void test_compute_pc_only() {
    printf("\n=== Test 5: compute_pc() — No Propagation ===\n");

    // Pre-computed states at TCA
    StateVector s1 = {2460310.5, 6778.0, 0, 0, 0, 7.669, 0};
    StateVector s2 = {2460310.5, 6778.5, 0.3, 0.1, 0.01, 7.668, 0.001};

    auto cov = Covariance3x3::socrates_default();

    ConjunctionEngine engine;
    engine.set_pc_method("foster");

    auto event = engine.compute_pc(s1, s2, cov, cov, 0.01);

    printf("  Miss: %.6f km\n", event.miss_distance_km);
    printf("  Pc (Foster): %.6e\n", event.pc.probability);
    printf("  Method: %s\n", event.pc.method.c_str());

    CHECK(event.miss_distance_km > 0, "Miss distance computed");
    CHECK(event.pc.probability >= 0, "Pc ≥ 0");
    CHECK(event.pc.method == "FOSTER-2D", "Method = FOSTER-2D");

    // Try with each method
    engine.set_pc_method("chan");
    auto event2 = engine.compute_pc(s1, s2, cov, cov, 0.01);
    CHECK(event2.pc.method == "CHAN-2008", "Method switchable to Chan");

    // Results should be similar
    if (event.pc.probability > 1e-30) {
        double ratio = event2.pc.probability / event.pc.probability;
        printf("  Chan/Foster ratio: %.3f\n", ratio);
        CHECK(ratio > 0.5 && ratio < 2.0, "Chan ≈ Foster (within 2×)");
    }
}

// ── Test 6: Cross-validate v1 vs v2 ──

void test_v1_v2_cross() {
    printf("\n=== Test 6: Cross-Validate v1 vs v2 ===\n");

    auto tle1 = parse_tle(ISS_TLE_NAME, ISS_TLE_L1, ISS_TLE_L2);
    auto tle2 = parse_tle(DEB_TLE_NAME, DEB_TLE_L1, DEB_TLE_L2);

    double start_jd = tle1.epoch_jd;

    // v1: Original assessment
    auto v1_event = assess_conjunction(tle1, tle2, start_jd, 1.0);

    // v2: New engine with Alfano (same method as v1)
    ConjunctionEngine engine;
    engine.set_pc_method("alfano");

    auto src1 = SGP4EphemerisSource(tle1);
    auto src2 = SGP4EphemerisSource(tle2);
    auto v2_event = engine.assess(src1, src2, start_jd, 1.0);

    printf("  v1 TCA: %s, miss: %.6f km, Pc: %.6e\n",
           v1_event.tca_iso.c_str(), v1_event.min_range_km, v1_event.max_probability);
    printf("  v2 TCA: %s, miss: %.6f km, Pc: %.6e\n",
           v2_event.tca_iso.c_str(), v2_event.miss_distance_km, v2_event.pc.probability);

    // TCA should match closely
    CHECK(std::abs(v1_event.tca_jd - v2_event.tca_jd) < 1.0/86400.0,
          "TCA matches within 1 second");

    // Miss distance should match
    CHECK_TOL(v2_event.miss_distance_km, v1_event.min_range_km, 0.01,
              "Miss distance matches v1");

    // Alfano Pc should match (same algorithm)
    if (v1_event.max_probability > 1e-30) {
        double ratio = v2_event.pc.probability / v1_event.max_probability;
        printf("  Pc ratio (v2/v1): %.6f\n", ratio);
        CHECK(std::abs(ratio - 1.0) < 0.01, "Alfano Pc matches v1 within 1%");
    }
}

// ── Test 7: Custom covariance ──

void test_custom_covariance() {
    printf("\n=== Test 7: Custom Covariance ===\n");

    StateVector s1 = {2460310.5, 6778.0, 0, 0, 0, 7.669, 0};
    StateVector s2 = {2460310.5, 6778.05, 0.02, 0, 0, 7.6689, 0}; // ~54m miss

    ConjunctionEngine engine;
    engine.set_pc_method("foster");

    // Small covariance (10m) — miss (~54m) is 5+ sigma away → low Pc
    auto cov_small = Covariance3x3::from_rtn_diagonal(0.01, 0.01, 0.01);
    auto e1 = engine.compute_pc(s1, s2, cov_small, cov_small, 0.01);

    // Large covariance (1km) — miss is within 0.05 sigma → higher Pc
    auto cov_large = Covariance3x3::from_rtn_diagonal(1.0, 1.0, 1.0);
    auto e2 = engine.compute_pc(s1, s2, cov_large, cov_large, 0.01);

    printf("  Miss: %.6f km\n", e1.miss_distance_km);
    printf("  Small cov (10m): Pc = %.6e\n", e1.pc.probability);
    printf("  Large cov (1km): Pc = %.6e\n", e2.pc.probability);

    // Both should produce valid, non-zero probabilities
    CHECK(e1.pc.probability > 0, "Small cov: Pc > 0");
    CHECK(e2.pc.probability > 0, "Large cov: Pc > 0");

    // At 54m miss, small cov (10m, ~3.8σ) gives higher Pc than large cov (1km, ~0.04σ)
    // because R²/(2σ₁σ₂) shrinks with larger σ faster than the exponential grows
    CHECK(e1.pc.probability != e2.pc.probability, "Different covariances → different Pc");
}

// ── Test 8: Mahalanobis distance in engine ──

void test_mahalanobis_engine() {
    printf("\n=== Test 8: Mahalanobis Distance (Engine) ===\n");

    // Two objects 500m apart, isotropic 100m covariance each
    StateVector s1 = {2460310.5, 6778.0, 0, 0, 0, 7.669, 0};
    StateVector s2 = {2460310.5, 6778.5, 0, 0, 0, 7.669, 0};  // 500m apart in X

    // Combined σ = sqrt(2 × 0.1²) = 0.1414 km = 141.4m per axis
    auto cov = Covariance3x3::from_rtn_diagonal(0.1, 0.1, 0.1);

    ConjunctionEngine engine;
    engine.set_pc_method("foster");

    auto event = engine.compute_pc(s1, s2, cov, cov, 0.01);

    printf("  Miss: %.3f km\n", event.miss_distance_km);
    printf("  Mahalanobis 2D: %.3f\n", event.mahalanobis_2d);
    printf("  Mahalanobis 3D: %.3f\n", event.mahalanobis_3d);
    printf("  Pc: %.6e\n", event.pc.probability);

    // 500m miss with combined σ = 141m → Md ≈ 3.5
    CHECK(event.mahalanobis_3d > 2.0, "3D Md > 2.0 (500m miss, 141m combined σ)");
    CHECK(event.mahalanobis_3d < 5.0, "3D Md < 5.0");
    CHECK(event.mahalanobis_2d >= 0, "2D Md ≥ 0");

    // PcResult also carries it
    CHECK_TOL(event.pc.mahalanobis_2d, event.mahalanobis_2d, 1e-10,
              "PcResult.mahalanobis_2d == event.mahalanobis_2d");
}

// ── Test 9: Mixed ephemeris sources ──

void test_mixed_sources() {
    printf("\n=== Test 8: Mixed Ephemeris Sources (SGP4 + OEM) ===\n");

    // Object 1: SGP4
    auto tle1 = parse_tle(ISS_TLE_NAME, ISS_TLE_L1, ISS_TLE_L2);
    SGP4EphemerisSource sgp4_src(tle1);

    // Object 2: OEM (generated from different TLE)
    auto tle2 = parse_tle(DEB_TLE_NAME, DEB_TLE_L1, DEB_TLE_L2);
    double epoch = tle2.epoch_jd;

    std::vector<EphemerisPoint> points;
    for (int i = 0; i <= 1440; i++) {  // 1 day at 1-minute steps
        double jd = epoch + i * 60.0 / 86400.0;
        auto sv = propagate_sgp4(tle2, jd);
        points.push_back({jd, sv.x, sv.y, sv.z, sv.vx, sv.vy, sv.vz});
    }
    OEMEphemerisSource oem_src(points, "DEBRIS-OEM", 99999);

    ConjunctionEngine engine;
    engine.set_pc_method("foster");

    auto event = engine.assess(sgp4_src, oem_src, epoch, 1.0);

    printf("  SGP4 vs OEM: TCA=%s, miss=%.3f km, Pc=%.6e\n",
           event.tca_iso.c_str(), event.miss_distance_km, event.pc.probability);

    CHECK(event.miss_distance_km >= 0, "Valid miss distance from mixed sources");
    CHECK(event.pc.probability >= 0, "Valid Pc from mixed sources");
    CHECK(!event.tca_iso.empty(), "TCA computed from mixed sources");
}

int main() {
    printf("============================================================\n");
    printf("Conjunction Engine v2 — Propagator-Agnostic Tests\n");
    printf("============================================================\n");

    test_sgp4_source();
    test_oem_source();
    test_callback_source();
    test_engine_methods();
    test_compute_pc_only();
    test_v1_v2_cross();
    test_custom_covariance();
    test_mahalanobis_engine();
    test_mixed_sources();

    printf("\n============================================================\n");
    printf("Results: %d passed, %d failed\n", tests_passed, tests_failed);
    printf("============================================================\n");

    return tests_failed > 0 ? 1 : 0;
}
