/**
 * Bulk SOCRATES Validation Test
 *
 * Tests conjunction assessment against hundreds/thousands of SOCRATES
 * conjunctions using the bulk CelesTrak GP JSON catalog.
 *
 * Usage:
 *   ./test_socrates_bulk <data_dir> [max_tests] [tolerance_mode]
 *
 * tolerance_mode:
 *   strict  — TCA < 10s, range < 100m, speed < 0.1 km/s
 *   normal  — TCA < 60s, range < 500m, speed < 0.5 km/s (default)
 *   relaxed — TCA < 300s, range < 5km, speed < 1.0 km/s
 */

#include "conjunction/conjunction_assessment.h"
#include "conjunction/gp_json.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <chrono>
#include <algorithm>
#include <iomanip>

using namespace conjunction;

// SOCRATES CSV record
struct SocratesRecord {
    int norad1, norad2;
    std::string name1, name2;
    double dse1, dse2;
    std::string tca;
    double min_range_km;
    double rel_speed_kms;
    double max_prob;
    double dilution_km;
};

// Parse SOCRATES CSV
static std::vector<SocratesRecord> parse_socrates_csv(const std::string& path, int max_records) {
    std::vector<SocratesRecord> records;
    std::ifstream f(path);
    if (!f.is_open()) return records;

    std::string line;
    std::getline(f, line); // Skip header

    while (std::getline(f, line) && (int)records.size() < max_records) {
        // Trim
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n'))
            line.pop_back();
        if (line.empty()) continue;

        // Parse CSV (handle quoted fields)
        std::vector<std::string> fields;
        std::string current;
        bool in_quotes = false;
        for (char c : line) {
            if (c == '"') { in_quotes = !in_quotes; continue; }
            if (c == ',' && !in_quotes) { fields.push_back(current); current.clear(); continue; }
            current += c;
        }
        fields.push_back(current);

        if (fields.size() < 11) continue;

        try {
            SocratesRecord rec;
            rec.norad1 = std::stoi(fields[0]);
            rec.name1 = fields[1];
            rec.dse1 = std::stod(fields[2]);
            rec.norad2 = std::stoi(fields[3]);
            rec.name2 = fields[4];
            rec.dse2 = std::stod(fields[5]);
            rec.tca = fields[6];
            // Normalize TCA format: "2026-03-14 11:08:50.281" → "2026-03-14T11:08:50.281Z"
            if (rec.tca.size() > 10 && rec.tca[10] == ' ') rec.tca[10] = 'T';
            if (rec.tca.back() != 'Z') rec.tca += 'Z';
            rec.min_range_km = std::stod(fields[7]);
            rec.rel_speed_kms = std::stod(fields[8]);
            rec.max_prob = std::stod(fields[9]);
            rec.dilution_km = std::stod(fields[10]);
            records.push_back(rec);
        } catch (...) {
            // Skip malformed lines
        }
    }

    return records;
}

// Load bulk GP catalog into a map by NORAD ID
static std::map<int, GPElement> load_gp_catalog(const std::string& path) {
    std::map<int, GPElement> catalog;
    std::ifstream f(path);
    if (!f.is_open()) return catalog;

    std::stringstream ss;
    ss << f.rdbuf();
    auto elements = parse_gp_json(ss.str());

    for (auto& gp : elements) {
        catalog[gp.norad_cat_id] = gp;
    }

    return catalog;
}

int main(int argc, char* argv[]) {
    std::string data_dir = argc > 1 ? argv[1] : "tests/data";
    int max_tests = argc > 2 ? std::atoi(argv[2]) : 100;
    std::string mode = argc > 3 ? argv[3] : "normal";

    // Tolerance settings
    double tca_tol_sec, range_tol_km, speed_tol_kms, prob_ratio_max;
    if (mode == "strict") {
        tca_tol_sec = 10.0; range_tol_km = 0.1; speed_tol_kms = 0.1; prob_ratio_max = 10.0;
    } else if (mode == "relaxed") {
        tca_tol_sec = 300.0; range_tol_km = 5.0; speed_tol_kms = 1.0; prob_ratio_max = 1000.0;
    } else {
        tca_tol_sec = 60.0; range_tol_km = 0.5; speed_tol_kms = 0.5; prob_ratio_max = 100.0;
    }

    std::cout << "=== Bulk SOCRATES Validation ===" << std::endl;
    std::cout << "Mode: " << mode << " (TCA<" << tca_tol_sec << "s, range<"
              << range_tol_km << "km, speed<" << speed_tol_kms << "km/s)" << std::endl;
    std::cout << "Max tests: " << max_tests << std::endl;

    // Load GP catalog (general, for fallback)
    std::cout << "\nLoading GP catalog..." << std::flush;
    auto catalog = load_gp_catalog(data_dir + "/celestrak_all_gp.json");
    std::cout << " " << catalog.size() << " objects" << std::endl;

    // SOCRATES GP directory (per-pair, exact epoch match)
    std::string socrates_gp_dir = data_dir + "/socrates_gp";

    // Load SOCRATES CSV (prefer current, fall back to maxprob)
    std::string csv_path = data_dir + "/socrates_current.csv";
    {
        std::ifstream test(csv_path);
        if (!test.is_open()) csv_path = data_dir + "/socrates_maxprob.csv";
    }
    std::cout << "Loading SOCRATES CSV: " << csv_path << "..." << std::flush;
    auto socrates = parse_socrates_csv(csv_path, max_tests * 2);
    std::cout << " " << socrates.size() << " records" << std::endl;

    // Each test uses assess_conjunction_near() with ±2 hour window around expected TCA

    int tested = 0, passed = 0, failed = 0, skipped = 0;
    int tca_pass = 0, range_pass = 0, speed_pass = 0, prob_pass = 0;

    // Track failure reasons
    std::map<std::string, int> failure_reasons;

    // Timing
    auto t_start = std::chrono::high_resolution_clock::now();
    double total_assess_ms = 0;

    std::cout << "\n--- Running " << max_tests << " tests ---" << std::endl;

    for (const auto& ref : socrates) {
        if (tested >= max_tests) break;

        // Skip formation-flight pairs (relative speed ≈ 0 means co-located constellation)
        // These have dozens of "conjunctions" per pair and are not actionable
        if (ref.rel_speed_kms < 0.01) {
            skipped++;
            failure_reasons["formation_flight"]++;
            continue;
        }

        // Alpha-5 encoding handles 6-digit NORAD IDs (100000-339999)

        // Try SOCRATES-epoch GP file first (exact match), fall back to catalog
        TLE tle1, tle2;
        bool loaded = false;

        // Try per-pair SOCRATES GP JSON
        std::string gp_path = socrates_gp_dir + "/gp_" +
            std::to_string(ref.norad1) + "," + std::to_string(ref.norad2) + ".json";
        {
            std::ifstream gp_file(gp_path);
            if (gp_file.is_open()) {
                std::stringstream ss; ss << gp_file.rdbuf();
                auto gps = parse_gp_json(ss.str());
                GPElement gp1_found, gp2_found;
                bool f1 = false, f2 = false;
                for (auto& gp : gps) {
                    if (gp.norad_cat_id == ref.norad1) { gp1_found = gp; f1 = true; }
                    if (gp.norad_cat_id == ref.norad2) { gp2_found = gp; f2 = true; }
                }
                if (f1 && f2) {
                    try {
                        tle1 = gp_to_tle(gp1_found);
                        tle2 = gp_to_tle(gp2_found);
                        loaded = true;
                    } catch (...) {}
                }
            }
        }

        // Fall back to general catalog
        if (!loaded) {
            auto it1 = catalog.find(ref.norad1);
            auto it2 = catalog.find(ref.norad2);
            if (it1 == catalog.end() || it2 == catalog.end()) {
                skipped++;
                continue;
            }
            try {
                tle1 = gp_to_tle(it1->second);
                tle2 = gp_to_tle(it2->second);
                loaded = true;
            } catch (...) {
                skipped++;
                failure_reasons["gp_to_tle_error"]++;
                continue;
            }
        }

        if (!loaded) { skipped++; continue; }

        tested++;

        // Run conjunction assessment near the expected TCA
        try {
            double ref_tca_jd_hint = iso_to_jd(ref.tca);
            auto t1 = std::chrono::high_resolution_clock::now();
            auto event = assess_conjunction_near(tle1, tle2, ref_tca_jd_hint, 0.75);
            auto t2 = std::chrono::high_resolution_clock::now();
            double assess_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
            total_assess_ms += assess_ms;

            // Compare
            double ref_tca_jd = iso_to_jd(ref.tca);
            double tca_diff = std::abs(event.tca_jd - ref_tca_jd) * 86400.0;
            double range_diff = std::abs(event.min_range_km - ref.min_range_km);
            double speed_diff = std::abs(event.rel_speed_kms - ref.rel_speed_kms);
            double prob_ratio = (ref.max_prob > 0 && event.max_probability > 0) ?
                std::max(event.max_probability / ref.max_prob,
                         ref.max_prob / event.max_probability) : 999;

            bool tca_ok = tca_diff < tca_tol_sec;
            bool range_ok = range_diff < range_tol_km;
            bool speed_ok = speed_diff < speed_tol_kms;
            bool prob_ok = prob_ratio < prob_ratio_max;

            if (tca_ok) tca_pass++;
            if (range_ok) range_pass++;
            if (speed_ok) speed_pass++;
            if (prob_ok) prob_pass++;

            bool all_ok = tca_ok && range_ok && speed_ok && prob_ok;

            if (all_ok) {
                passed++;
            } else {
                failed++;

                // Categorize failure
                if (!tca_ok) failure_reasons["tca_mismatch"]++;
                if (!range_ok) failure_reasons["range_mismatch"]++;
                if (!speed_ok) failure_reasons["speed_mismatch"]++;
                if (!prob_ok) failure_reasons["prob_mismatch"]++;

                // Print failures (first 20 + every 50th)
                if (failed <= 20 || failed % 50 == 0) {
                    std::cout << "  FAIL #" << tested << " [" << ref.norad1 << " vs " << ref.norad2 << "]"
                              << " TCA_Δ=" << std::fixed << std::setprecision(1) << tca_diff << "s"
                              << " range=" << std::setprecision(3) << event.min_range_km
                              << "/" << ref.min_range_km << "km"
                              << " speed=" << std::setprecision(3) << event.rel_speed_kms
                              << "/" << ref.rel_speed_kms << "km/s"
                              << " prob=" << std::scientific << std::setprecision(2) << event.max_probability
                              << "/" << ref.max_prob
                              << " (" << std::fixed << std::setprecision(0) << assess_ms << "ms)"
                              << std::endl;
                }
            }

            // Progress every 50 tests
            if (tested % 50 == 0) {
                std::cout << "  [" << tested << "/" << max_tests << "] "
                          << passed << " pass, " << failed << " fail, " << skipped << " skip"
                          << " (avg " << std::fixed << std::setprecision(0)
                          << total_assess_ms / tested << "ms/test)" << std::endl;
            }

        } catch (const std::exception& e) {
            failed++;
            failure_reasons["exception"]++;
            if (failed <= 20) {
                std::cout << "  ERROR #" << tested << " [" << ref.norad1 << " vs " << ref.norad2 << "] "
                          << e.what() << std::endl;
            }
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    // Summary
    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Tested:  " << tested << std::endl;
    std::cout << "Passed:  " << passed << " (" << std::fixed << std::setprecision(1)
              << (tested > 0 ? 100.0 * passed / tested : 0) << "%)" << std::endl;
    std::cout << "Failed:  " << failed << " (" << (tested > 0 ? 100.0 * failed / tested : 0) << "%)" << std::endl;
    std::cout << "Skipped: " << skipped << " (missing GP data)" << std::endl;

    std::cout << "\nPer-metric pass rate:" << std::endl;
    std::cout << "  TCA:   " << tca_pass << "/" << tested << " (" << (tested > 0 ? 100.0 * tca_pass / tested : 0) << "%)" << std::endl;
    std::cout << "  Range: " << range_pass << "/" << tested << " (" << (tested > 0 ? 100.0 * range_pass / tested : 0) << "%)" << std::endl;
    std::cout << "  Speed: " << speed_pass << "/" << tested << " (" << (tested > 0 ? 100.0 * speed_pass / tested : 0) << "%)" << std::endl;
    std::cout << "  Prob:  " << prob_pass << "/" << tested << " (" << (tested > 0 ? 100.0 * prob_pass / tested : 0) << "%)" << std::endl;

    if (!failure_reasons.empty()) {
        std::cout << "\nFailure reasons:" << std::endl;
        for (auto& [reason, count] : failure_reasons) {
            std::cout << "  " << reason << ": " << count << std::endl;
        }
    }

    std::cout << "\nTiming:" << std::endl;
    std::cout << "  Total: " << std::fixed << std::setprecision(1) << total_ms << " ms" << std::endl;
    std::cout << "  Avg per test: " << (tested > 0 ? total_ms / tested : 0) << " ms" << std::endl;
    std::cout << "  Avg assess: " << (tested > 0 ? total_assess_ms / tested : 0) << " ms" << std::endl;

    return (failed > 0) ? 1 : 0;
}
