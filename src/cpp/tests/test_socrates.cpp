/**
 * Test: Validate conjunction assessment against SOCRATES data from CelesTrak
 *
 * Loads TLE pairs from SOCRATES GP data files and compares computed
 * TCA, miss distance, relative speed, and max probability against
 * SOCRATES reference values.
 */

#include "conjunction/conjunction_assessment.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>

using namespace conjunction;

// SOCRATES reference data (scraped 2026-03-10)
struct SocratesRef {
    const char* gp_file;
    int norad1, norad2;
    const char* tca;           // ISO 8601
    double min_range_km;
    double rel_speed_kms;
    double max_prob;
    double dilution_km;
    double dse1, dse2;
};

static const SocratesRef SOCRATES_DATA[] = {
    {"gp_61721,67298.txt", 61721, 67298, "2026-03-14T11:08:50.281Z",
     0.010, 9.757, 1.000, 0.000, 5.075, 5.733},
    {"gp_47935,49179.txt", 47935, 49179, "2026-03-12T04:44:40.733Z",
     0.014, 9.173, 0.3343, 0.004, 2.681, 2.684},
    {"gp_48282,58288.txt", 48282, 58288, "2026-03-16T09:03:22.330Z",
     0.014, 12.349, 0.3310, 0.005, 7.136, 6.872},
    {"gp_46054,61781.txt", 46054, 61781, "2026-03-15T23:13:06.039Z",
     0.013, 14.013, 0.3235, 0.005, 6.467, 6.417},
    {"gp_58184,54940.txt", 58184, 54940, "2026-03-15T01:07:39.931Z",
     0.029, 8.865, 0.1153, 0.008, 6.380, 5.811},
    {"gp_56844,64555.txt", 56844, 64555, "2026-03-14T01:29:50.899Z",
     0.030, 12.230, 0.09097, 0.011, 4.562, 4.746},
    {"gp_46335,59761.txt", 46335, 59761, "2026-03-14T19:57:29.280Z",
     0.061, 11.147, 0.07876, 0.020, 5.332, 5.278},
    {"gp_60150,43771.txt", 60150, 43771, "2026-03-10T23:11:36.576Z",
     0.031, 14.757, 0.05748, 0.017, 1.567, 2.321},
    {"gp_57446,59806.txt", 57446, 59806, "2026-03-11T18:37:48.752Z",
     0.046, 11.393, 0.04243, 0.015, 2.276, 2.269},
    {"gp_64428,20608.txt", 64428, 20608, "2026-03-15T04:50:51.304Z",
     0.042, 14.912, 0.03955, 0.026, 5.745, 5.707},
};

static std::string read_file(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) throw std::runtime_error("Cannot open: " + path);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

int main(int argc, char* argv[]) {
    std::string data_dir = "tests/data";
    if (argc > 1) data_dir = argv[1];

    // SOCRATES computation start
    double start_jd = iso_to_jd("2026-03-09T18:00:00Z");

    int passed = 0, failed = 0, errors = 0;
    const int N = sizeof(SOCRATES_DATA) / sizeof(SOCRATES_DATA[0]);

    std::cout << "=== SOCRATES Validation Test ===" << std::endl;
    std::cout << "Testing " << N << " conjunction events" << std::endl;
    std::cout << std::endl;

    for (int i = 0; i < N; i++) {
        const auto& ref = SOCRATES_DATA[i];
        std::cout << "--- Conjunction " << (i + 1) << ": "
                  << ref.norad1 << " vs " << ref.norad2 << " ---" << std::endl;

        try {
            // Load GP data
            std::string gp = read_file(data_dir + "/" + ref.gp_file);
            auto tles = parse_tle_file(gp);

            if (tles.size() < 2) {
                std::cerr << "  ERROR: Expected 2 TLEs, got " << tles.size() << std::endl;
                errors++;
                continue;
            }

            // Find the right TLEs by NORAD ID
            TLE tle1, tle2;
            bool found1 = false, found2 = false;
            for (const auto& t : tles) {
                if (t.norad_cat_id == ref.norad1) { tle1 = t; found1 = true; }
                if (t.norad_cat_id == ref.norad2) { tle2 = t; found2 = true; }
            }

            if (!found1 || !found2) {
                std::cerr << "  ERROR: Could not find TLEs for NORAD IDs" << std::endl;
                errors++;
                continue;
            }

            // Run conjunction assessment
            auto event = assess_conjunction(tle1, tle2, start_jd, 7.0);

            // Compare results
            double ref_tca_jd = iso_to_jd(ref.tca);
            double tca_diff_sec = std::abs(event.tca_jd - ref_tca_jd) * 86400.0;
            double range_diff_km = std::abs(event.min_range_km - ref.min_range_km);
            double speed_diff = std::abs(event.rel_speed_kms - ref.rel_speed_kms);

            // Tolerances
            bool tca_ok = tca_diff_sec < 60.0;       // < 1 minute
            bool range_ok = range_diff_km < 0.5;      // < 500m
            bool speed_ok = speed_diff < 0.5;          // < 0.5 km/s

            // Max probability comparison (order of magnitude)
            double prob_ratio = (ref.max_prob > 0 && event.max_probability > 0) ?
                std::max(event.max_probability / ref.max_prob,
                         ref.max_prob / event.max_probability) : 999;
            bool prob_ok = prob_ratio < 100.0;  // Within 2 orders of magnitude

            std::cout << "  TCA: " << event.tca_iso
                      << " (SOCRATES: " << ref.tca << ")"
                      << " Δ=" << std::fixed << std::setprecision(1) << tca_diff_sec << "s"
                      << (tca_ok ? " ✓" : " ✗") << std::endl;

            std::cout << "  Range: " << std::setprecision(3) << event.min_range_km
                      << " km (SOCRATES: " << ref.min_range_km << " km)"
                      << " Δ=" << range_diff_km << " km"
                      << (range_ok ? " ✓" : " ✗") << std::endl;

            std::cout << "  Speed: " << std::setprecision(3) << event.rel_speed_kms
                      << " km/s (SOCRATES: " << ref.rel_speed_kms << " km/s)"
                      << " Δ=" << speed_diff << " km/s"
                      << (speed_ok ? " ✓" : " ✗") << std::endl;

            std::cout << "  MaxProb: " << std::scientific << std::setprecision(3) << event.max_probability
                      << " (SOCRATES: " << ref.max_prob << ")"
                      << " ratio=" << std::fixed << std::setprecision(1) << prob_ratio
                      << (prob_ok ? " ✓" : " ✗") << std::endl;

            std::cout << "  DSE1: " << std::setprecision(3) << event.dse1
                      << " (SOCRATES: " << ref.dse1 << ")" << std::endl;
            std::cout << "  DSE2: " << std::setprecision(3) << event.dse2
                      << " (SOCRATES: " << ref.dse2 << ")" << std::endl;

            if (tca_ok && range_ok && speed_ok && prob_ok) {
                std::cout << "  RESULT: PASS ✅" << std::endl;
                passed++;
            } else {
                std::cout << "  RESULT: FAIL ❌" << std::endl;
                failed++;
            }

        } catch (const std::exception& e) {
            std::cerr << "  ERROR: " << e.what() << std::endl;
            errors++;
        }

        std::cout << std::endl;
    }

    std::cout << "=== Summary ===" << std::endl;
    std::cout << "Passed: " << passed << "/" << N << std::endl;
    std::cout << "Failed: " << failed << "/" << N << std::endl;
    std::cout << "Errors: " << errors << "/" << N << std::endl;

    return (failed + errors > 0) ? 1 : 0;
}
