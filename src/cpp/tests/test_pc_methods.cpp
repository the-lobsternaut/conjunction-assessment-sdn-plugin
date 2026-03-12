/**
 * Collision Probability Methods — Cross-Validation Tests
 *
 * Strategy:
 *   1. Known analytical cases (head-on, zero miss, large miss)
 *   2. Cross-method consistency (all methods should agree within bounds)
 *   3. Alfano always ≥ all other methods (it's the maximum)
 *   4. Convergence: high-resolution methods match each other
 *   5. Edge cases: degenerate covariance, circular covariance
 */

#include "conjunction/pc_method.h"
#include <cstdio>
#include <cmath>
#include <memory>
#include <vector>

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

// ── Test 1: All methods on a standard case ──

void test_standard_case() {
    printf("\n=== Test 1: Standard Case (1 km miss, 10m radius, 100m/300m/100m cov) ===\n");

    BPlaneGeometry bp;
    bp.xi = 0.7;    // 0.7 km miss in X
    bp.zeta = 0.7;  // 0.7 km miss in Z → ~1 km total miss
    bp.sigma_xx = 0.1 * 0.1 + 0.1 * 0.1;  // combined: 2 × (100m)² in km²
    bp.sigma_xz = 0;
    bp.sigma_zz = 0.3 * 0.3 + 0.3 * 0.3;  // combined: 2 × (300m)² in km²
    bp.combined_radius = 0.01;  // 10m

    printf("  Miss: %.3f km, Radius: %.3f km\n", bp.miss_distance(), bp.combined_radius);
    printf("  Cov eigenvalues: σ_xx=%.4f, σ_zz=%.4f km²\n", bp.sigma_xx, bp.sigma_zz);

    auto alfano = std::make_unique<AlfanoMaxPc>();
    auto foster = std::make_unique<Foster2D>();
    auto patera = std::make_unique<Patera2001>(256);
    auto chan = std::make_unique<Chan2008>();
    auto alfriend = std::make_unique<Alfriend2D>(64, 128);

    auto r_alf = alfano->compute(bp);
    auto r_fos = foster->compute(bp);
    auto r_pat = patera->compute(bp);
    auto r_cha = chan->compute(bp);
    auto r_afr = alfriend->compute(bp);

    printf("  Alfano:   %.6e\n", r_alf.probability);
    printf("  Foster:   %.6e\n", r_fos.probability);
    printf("  Patera:   %.6e\n", r_pat.probability);
    printf("  Chan:     %.6e\n", r_cha.probability);
    printf("  Alfriend: %.6e\n", r_afr.probability);

    // Alfano must be >= all others (it's the maximum probability)
    CHECK(r_alf.probability >= r_fos.probability * 0.99, "Alfano ≥ Foster");
    CHECK(r_alf.probability >= r_pat.probability * 0.99, "Alfano ≥ Patera");
    CHECK(r_alf.probability >= r_cha.probability * 0.99, "Alfano ≥ Chan");
    CHECK(r_alf.probability >= r_afr.probability * 0.99, "Alfano ≥ Alfriend");

    // Chan and Alfriend should agree closely (both do 2D numerical integration)
    if (r_cha.probability > 1e-30 && r_afr.probability > 1e-30) {
        double ratio_ca = r_cha.probability / r_afr.probability;
        printf("  Chan/Alfriend ratio: %.4f\n", ratio_ca);
        CHECK(ratio_ca > 0.9 && ratio_ca < 1.1, "Chan ≈ Alfriend within 10%");
    }

    // All should be non-zero
    CHECK(r_fos.probability > 0, "Foster > 0");
    CHECK(r_pat.probability > 0, "Patera > 0");
    CHECK(r_cha.probability > 0, "Chan > 0");
    CHECK(r_afr.probability > 0, "Alfriend > 0");
}

// ── Test 2: Zero miss distance ──

void test_zero_miss() {
    printf("\n=== Test 2: Zero Miss Distance ===\n");

    BPlaneGeometry bp;
    bp.xi = 0; bp.zeta = 0;
    bp.sigma_xx = 0.01; bp.sigma_zz = 0.01;
    bp.combined_radius = 0.01;

    auto foster = std::make_unique<Foster2D>();
    auto chan = std::make_unique<Chan2008>();

    auto r_fos = foster->compute(bp);
    auto r_cha = chan->compute(bp);

    printf("  Foster Pc: %.6e\n", r_fos.probability);
    printf("  Chan Pc:   %.6e\n", r_cha.probability);

    // At zero miss, Pc should be relatively high
    CHECK(r_fos.probability > 1e-6, "Foster: non-negligible Pc at zero miss");
    CHECK(r_cha.probability > 1e-6, "Chan: non-negligible Pc at zero miss");
}

// ── Test 3: Large miss distance ──

void test_large_miss() {
    printf("\n=== Test 3: Large Miss Distance (100 km) ===\n");

    BPlaneGeometry bp;
    bp.xi = 70.71; bp.zeta = 70.71;  // ~100 km miss
    bp.sigma_xx = 0.01; bp.sigma_zz = 0.01;
    bp.combined_radius = 0.01;

    auto foster = std::make_unique<Foster2D>();
    auto alfano = std::make_unique<AlfanoMaxPc>();

    auto r_fos = foster->compute(bp);
    auto r_alf = alfano->compute(bp);

    printf("  Foster Pc: %.6e\n", r_fos.probability);
    printf("  Alfano Pc: %.6e\n", r_alf.probability);

    // Both should be extremely small
    CHECK(r_fos.probability < 1e-20, "Foster: negligible at 100 km miss");
    CHECK(r_alf.probability < 1e-4, "Alfano: small at 100 km miss");
}

// ── Test 4: Circular covariance (Foster exact) ──

void test_circular_cov() {
    printf("\n=== Test 4: Circular Covariance (isotropic, Foster exact) ===\n");

    BPlaneGeometry bp;
    bp.xi = 0.5; bp.zeta = 0;
    bp.sigma_xx = 0.04; bp.sigma_zz = 0.04;  // σ = 200m isotropic
    bp.combined_radius = 0.01;

    auto foster = std::make_unique<Foster2D>();
    auto patera = std::make_unique<Patera2001>(512);
    auto chan = std::make_unique<Chan2008>();
    auto alfriend = std::make_unique<Alfriend2D>(64, 256);

    auto r_fos = foster->compute(bp);
    auto r_pat = patera->compute(bp);
    auto r_cha = chan->compute(bp);
    auto r_afr = alfriend->compute(bp);

    printf("  Foster:   %.10e\n", r_fos.probability);
    printf("  Patera:   %.10e\n", r_pat.probability);
    printf("  Chan:     %.10e\n", r_cha.probability);
    printf("  Alfriend: %.10e\n", r_afr.probability);

    // For circular covariance, all methods should agree closely
    // Use Chan/Alfriend (full integration) as reference
    double ref = r_afr.probability;
    if (ref > 1e-30) {
        CHECK(std::abs(r_fos.probability - ref) / ref < 0.05, "Foster within 5% of Alfriend (circular)");
        CHECK(std::abs(r_cha.probability - ref) / ref < 0.05, "Chan within 5% of Alfriend (circular)");
        CHECK(std::abs(r_pat.probability - ref) / ref < 0.1, "Patera within 10% of Alfriend (circular)");
    }
}

// ── Test 5: Highly elliptical covariance ──

void test_elliptical_cov() {
    printf("\n=== Test 5: Highly Elliptical Covariance (10:1 ratio) ===\n");

    BPlaneGeometry bp;
    bp.xi = 0.3; bp.zeta = 0.1;
    bp.sigma_xx = 0.001;   // σ = 31.6 m (tight)
    bp.sigma_zz = 0.1;     // σ = 316 m (wide) → 10:1 ratio
    bp.combined_radius = 0.01;

    auto foster = std::make_unique<Foster2D>();
    auto patera = std::make_unique<Patera2001>(512);
    auto alfriend = std::make_unique<Alfriend2D>(64, 256);

    auto r_fos = foster->compute(bp);
    auto r_pat = patera->compute(bp);
    auto r_afr = alfriend->compute(bp);

    printf("  Foster:   %.10e\n", r_fos.probability);
    printf("  Patera:   %.10e\n", r_pat.probability);
    printf("  Alfriend: %.10e\n", r_afr.probability);

    // All numerical methods should agree within an order of magnitude
    // for highly elliptical cases (exact values differ due to quadrature resolution)
    if (r_afr.probability > 1e-30 && r_pat.probability > 1e-30) {
        double ratio = r_pat.probability / r_afr.probability;
        printf("  Patera/Alfriend ratio: %.3f\n", ratio);
        CHECK(ratio > 0.1 && ratio < 10.0, "Patera ≈ Alfriend within 1 order of magnitude (elliptical)");
    }
}

// ── Test 6: Factory ──

void test_factory() {
    printf("\n=== Test 6: Method Factory ===\n");

    auto m1 = create_pc_method("alfano");
    CHECK(m1->name() == "ALFANO-MAXPROB", "Factory: alfano");

    auto m2 = create_pc_method("foster");
    CHECK(m2->name() == "FOSTER-2D", "Factory: foster");

    auto m3 = create_pc_method("patera");
    CHECK(m3->name() == "PATERA-2001", "Factory: patera");

    auto m4 = create_pc_method("chan");
    CHECK(m4->name() == "CHAN-2008", "Factory: chan");

    auto m5 = create_pc_method("alfriend");
    CHECK(m5->name() == "ALFRIEND-2D", "Factory: alfriend");

    auto m6 = create_pc_method("unknown");
    CHECK(m6->name() == "FOSTER-2D", "Factory: unknown → default Foster");
}

// ── Test 8: Mahalanobis Distance ──

void test_mahalanobis() {
    printf("\n=== Test 8: Mahalanobis Distance ===\n");

    // Case 1: Miss at 1 sigma in isotropic covariance
    BPlaneGeometry bp1;
    bp1.xi = 0.1; bp1.zeta = 0;
    bp1.sigma_xx = 0.01;  // σ = 100m
    bp1.sigma_xz = 0;
    bp1.sigma_zz = 0.01;  // σ = 100m
    bp1.combined_radius = 0.01;

    double md1 = bp1.mahalanobis_distance();
    printf("  100m miss, 100m σ isotropic: Md = %.3f\n", md1);
    CHECK_TOL(md1, 1.0, 0.01, "Md = 1.0 at 1σ miss (isotropic)");

    // Case 2: Miss at 3 sigma
    BPlaneGeometry bp2;
    bp2.xi = 0.3; bp2.zeta = 0;
    bp2.sigma_xx = 0.01; bp2.sigma_zz = 0.01;
    double md2 = bp2.mahalanobis_distance();
    printf("  300m miss, 100m σ: Md = %.3f\n", md2);
    CHECK_TOL(md2, 3.0, 0.01, "Md = 3.0 at 3σ miss");

    // Case 3: Elliptical covariance — miss along minor axis
    BPlaneGeometry bp3;
    bp3.xi = 0.1; bp3.zeta = 0;
    bp3.sigma_xx = 0.01;   // σ_x = 100m (tight)
    bp3.sigma_xz = 0;
    bp3.sigma_zz = 1.0;    // σ_z = 1km (wide)
    double md3 = bp3.mahalanobis_distance();
    printf("  100m miss along tight axis (σ=100m, σ=1km): Md = %.3f\n", md3);
    CHECK_TOL(md3, 1.0, 0.01, "Md = 1.0 along minor axis");

    // Case 4: Same miss along wide axis
    BPlaneGeometry bp4;
    bp4.xi = 0; bp4.zeta = 0.1;  // miss along wide axis
    bp4.sigma_xx = 0.01;
    bp4.sigma_xz = 0;
    bp4.sigma_zz = 1.0;
    double md4 = bp4.mahalanobis_distance();
    printf("  100m miss along wide axis (σ=1km): Md = %.3f\n", md4);
    CHECK(md4 < 0.15, "Md < 0.15 (miss along wide axis, well within 1σ)");

    // Case 5: Zero miss
    BPlaneGeometry bp5;
    bp5.xi = 0; bp5.zeta = 0;
    bp5.sigma_xx = 0.01; bp5.sigma_zz = 0.01;
    double md5 = bp5.mahalanobis_distance();
    CHECK_TOL(md5, 0.0, 1e-10, "Md = 0 at zero miss");

    // Case 6: Mahalanobis stored in PcResult
    auto foster = std::make_unique<Foster2D>();
    auto r = foster->compute(bp1);
    printf("  PcResult.mahalanobis_2d = %.3f\n", r.mahalanobis_2d);
    CHECK_TOL(r.mahalanobis_2d, 1.0, 0.01, "PcResult carries Mahalanobis distance");
}

// ── Test 7: Overlapping hard bodies ──

void test_overlapping() {
    printf("\n=== Test 7: Overlapping Hard Bodies ===\n");

    BPlaneGeometry bp;
    bp.xi = 0.005; bp.zeta = 0;  // 5m miss
    bp.sigma_xx = 0.01; bp.sigma_zz = 0.01;
    bp.combined_radius = 0.01;  // 10m radius > 5m miss

    auto alfano = std::make_unique<AlfanoMaxPc>();
    auto r = alfano->compute(bp);

    CHECK(r.probability == 1.0, "Alfano: Pc = 1.0 when radius > miss");
}

int main() {
    printf("============================================================\n");
    printf("Collision Probability Methods — Cross-Validation\n");
    printf("============================================================\n");

    test_standard_case();
    test_zero_miss();
    test_large_miss();
    test_circular_cov();
    test_elliptical_cov();
    test_factory();
    test_overlapping();
    test_mahalanobis();

    printf("\n============================================================\n");
    printf("Results: %d passed, %d failed\n", tests_passed, tests_failed);
    printf("============================================================\n");

    return tests_failed > 0 ? 1 : 0;
}
