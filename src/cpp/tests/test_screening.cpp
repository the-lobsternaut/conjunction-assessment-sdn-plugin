/**
 * Test: KD-Tree Screening Engine Benchmark
 *
 * Tests the high-performance screening pipeline:
 *   - KD-tree build + range query correctness
 *   - Perigee/apogee prefilter
 *   - Multi-threaded screening
 *   - Performance comparison: KD-tree vs brute force
 */

#include "conjunction/screening.h"
#include "conjunction/gp_json.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cmath>
#include <random>

using namespace conjunction;

static std::string read_file(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) throw std::runtime_error("Cannot open: " + path);
    std::stringstream ss; ss << f.rdbuf();
    return ss.str();
}

void test_kdtree_correctness() {
    std::cout << "\n=== KD-Tree Correctness ===" << std::endl;

    // Create a set of known points
    std::vector<KDPoint> points = {
        {100.0, 200.0, 300.0, 0, 0},
        {101.0, 200.0, 300.0, 1, 0},  // 1 km from point 0
        {200.0, 300.0, 400.0, 2, 0},  // Far from 0 and 1
        {100.5, 200.5, 300.5, 3, 0},  // ~0.87 km from point 0
    };

    KDTree tree;
    tree.build(points);

    // Query for neighbors of point 0 within 2 km
    std::vector<uint32_t> results;
    tree.range_query(points[0], 2.0, results);

    std::cout << "  Points within 2 km of (100,200,300): " << results.size() << std::endl;
    bool found_self = false, found_1 = false, found_3 = false, found_2 = false;
    for (auto idx : results) {
        auto& p = points[idx];
        if (p.obj_index == 0) found_self = true;
        if (p.obj_index == 1) found_1 = true;
        if (p.obj_index == 2) found_2 = true;
        if (p.obj_index == 3) found_3 = true;
    }

    std::cout << "  Found self: " << (found_self ? "✓" : "✗") << std::endl;
    std::cout << "  Found point 1 (1 km): " << (found_1 ? "✓" : "✗") << std::endl;
    std::cout << "  Found point 3 (0.87 km): " << (found_3 ? "✓" : "✗") << std::endl;
    std::cout << "  Excluded point 2 (far): " << (!found_2 ? "✓" : "✗") << std::endl;
}

void test_kdtree_performance() {
    std::cout << "\n=== KD-Tree Performance ===" << std::endl;

    // Generate N random LEO positions
    std::mt19937 rng(42);
    std::uniform_real_distribution<> angle_dist(0, 6.283185);
    std::normal_distribution<> r_dist(6778.0, 100.0); // ~400km altitude

    std::vector<int> sizes = {100, 500, 1000, 5000};

    for (int N : sizes) {
        std::vector<KDPoint> points(N);
        for (int i = 0; i < N; i++) {
            double r = r_dist(rng);
            double theta = angle_dist(rng);
            double phi = angle_dist(rng);
            points[i] = {
                r * std::sin(theta) * std::cos(phi),
                r * std::sin(theta) * std::sin(phi),
                r * std::cos(theta),
                static_cast<uint32_t>(i), 0
            };
        }

        // Time KD-tree build
        auto t1 = std::chrono::high_resolution_clock::now();
        KDTree tree;
        tree.build(points);
        auto t2 = std::chrono::high_resolution_clock::now();
        double build_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();

        // Time KD-tree find_close_pairs (5 km threshold)
        auto t3 = std::chrono::high_resolution_clock::now();
        auto pairs = tree.find_close_pairs(points, 5.0);
        auto t4 = std::chrono::high_resolution_clock::now();
        double query_ms = std::chrono::duration<double, std::milli>(t4 - t3).count();

        // Time brute force for comparison
        auto t5 = std::chrono::high_resolution_clock::now();
        int brute_count = 0;
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                double dx = points[i].x - points[j].x;
                double dy = points[i].y - points[j].y;
                double dz = points[i].z - points[j].z;
                if (dx*dx + dy*dy + dz*dz <= 25.0) brute_count++;
            }
        }
        auto t6 = std::chrono::high_resolution_clock::now();
        double brute_ms = std::chrono::duration<double, std::milli>(t6 - t5).count();

        std::cout << "  N=" << N
                  << " | build=" << std::fixed << std::setprecision(1) << build_ms << "ms"
                  << " | query=" << query_ms << "ms"
                  << " | brute=" << brute_ms << "ms"
                  << " | speedup=" << std::setprecision(1) << (brute_ms / std::max(0.001, query_ms)) << "x"
                  << " | pairs=" << pairs.size()
                  << " (brute=" << brute_count << ")"
                  << std::endl;
    }
}

void test_screening_with_socrates(const std::string& data_dir) {
    std::cout << "\n=== Full Screening Pipeline (SOCRATES top 10) ===" << std::endl;

    // Load SOCRATES JSON fixture
    std::string fixture_path = data_dir + "/socrates_top20_maxProb.json";
    std::ifstream f(fixture_path);
    if (!f.is_open()) {
        std::cout << "  Skipped: " << fixture_path << " not found" << std::endl;
        return;
    }

    // Load all GP data as JSON
    std::vector<GPElement> catalog;
    for (int i = 1; i <= 10; i++) {
        // Read the first 10 GP files
        std::string gp_file = data_dir + "/gp_";
        // We need to know the NORAD IDs — read from the fixture
        // For now, test with a simple pair
    }

    // Simple test: load one GP JSON pair and screen
    const char* json = R"([
        {"OBJECT_NAME":"STARLINK-32469","OBJECT_ID":"2024-196Q",
         "EPOCH":"2026-03-09T09:20:49.941888","MEAN_MOTION":15.30192485,
         "ECCENTRICITY":0.0001054,"INCLINATION":53.1601,
         "RA_OF_ASC_NODE":102.5175,"ARG_OF_PERICENTER":81.5233,
         "MEAN_ANOMALY":278.5887,"EPHEMERIS_TYPE":0,
         "CLASSIFICATION_TYPE":"U","NORAD_CAT_ID":61721,
         "ELEMENT_SET_NO":999,"REV_AT_EPOCH":7711,
         "BSTAR":3.1791e-5,"MEAN_MOTION_DOT":5.82e-6,"MEAN_MOTION_DDOT":0},
        {"OBJECT_NAME":"OBJECT BC","OBJECT_ID":"2025-313BC",
         "EPOCH":"2026-03-08T17:33:30.106080","MEAN_MOTION":15.29053282,
         "ECCENTRICITY":0.00085081,"INCLINATION":97.3987,
         "RA_OF_ASC_NODE":144.3551,"ARG_OF_PERICENTER":331.493,
         "MEAN_ANOMALY":28.5844,"EPHEMERIS_TYPE":0,
         "CLASSIFICATION_TYPE":"U","NORAD_CAT_ID":67298,
         "ELEMENT_SET_NO":999,"REV_AT_EPOCH":782,
         "BSTAR":0.0039787182,"MEAN_MOTION_DOT":0.00113849,"MEAN_MOTION_DDOT":0}
    ])";

    auto elements = parse_gp_json(json);

    ScreeningConfig config;
    config.start_jd = iso_to_jd("2026-03-09T18:00:00Z");
    config.duration_days = 7.0;
    config.threshold_km = 5.0;
    config.num_threads = 4;
    config.use_kdtree = false; // Only 2 objects, no point
    config.coarse_step_sec = 60.0;

    ConjunctionScreener screener(config);
    auto events = screener.screen(elements,
        [](double frac, const std::string& status) {
            std::cout << "  [" << std::fixed << std::setprecision(0)
                      << (frac * 100) << "%] " << status << std::endl;
        });

    auto& stats = screener.stats();
    std::cout << "\n  Stats:" << std::endl;
    std::cout << "    Objects: " << stats.total_objects << std::endl;
    std::cout << "    Pairs screened: " << stats.pairs_screened << std::endl;
    std::cout << "    Prefiltered: " << stats.pairs_prefiltered << std::endl;
    std::cout << "    KD-tree candidates: " << stats.kdtree_candidates << std::endl;
    std::cout << "    TCA refined: " << stats.tca_refined << std::endl;
    std::cout << "    Conjunctions found: " << stats.conjunctions_found << std::endl;
    std::cout << "    Propagations: " << stats.propagations << std::endl;
    std::cout << "    Elapsed: " << std::fixed << std::setprecision(1)
              << stats.elapsed_ms << " ms" << std::endl;

    if (!events.empty()) {
        std::cout << "\n  Top conjunction:" << std::endl;
        std::cout << "    TCA: " << events[0].tca_iso << std::endl;
        std::cout << "    Range: " << std::setprecision(3)
                  << events[0].min_range_km << " km" << std::endl;
        std::cout << "    Speed: " << events[0].rel_speed_kms << " km/s" << std::endl;
        std::cout << "    MaxProb: " << std::scientific
                  << events[0].max_probability << std::endl;
    }
}

int main(int argc, char* argv[]) {
    std::string data_dir = "tests/data";
    if (argc > 1) data_dir = argv[1];

    test_kdtree_correctness();
    test_kdtree_performance();
    test_screening_with_socrates(data_dir);

    std::cout << "\n=== Done ===" << std::endl;
    return 0;
}
