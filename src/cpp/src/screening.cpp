/**
 * High-Performance Conjunction Screening Engine
 *
 * Pipeline:
 *   1. Perigee/apogee prefilter (eliminate impossible pairs)
 *   2. Convert GP → TLE for SGP4 propagation
 *   3. Parallel time-step propagation with KD-tree spatial indexing
 *   4. Dynamic windowing (adaptive step based on closing rate)
 *   5. Fine TCA refinement for candidates
 *   6. Full conjunction assessment for confirmed events
 *
 * Threading: pthreads with per-thread work queues
 * Memory: shared read-only TLE array, per-thread candidate lists
 */

#include "conjunction/screening.h"
#include <thread>
#include <mutex>
#include <chrono>
#include <cmath>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <map>
#include <iostream>

namespace conjunction {

// ============================================================================
// Perigee/Apogee Prefilter
// ============================================================================

bool altitude_overlap(const GPElement& gp1, const GPElement& gp2, double margin_km) {
    // Objects can only conjunct if their altitude ranges overlap
    // Perigee1-margin .. Apogee1+margin overlaps with Perigee2-margin .. Apogee2+margin
    double lo1 = gp1.perigee_km - margin_km;
    double hi1 = gp1.apogee_km + margin_km;
    double lo2 = gp2.perigee_km - margin_km;
    double hi2 = gp2.apogee_km + margin_km;

    return lo1 <= hi2 && lo2 <= hi1;
}

std::vector<std::pair<uint32_t, uint32_t>> ConjunctionScreener::prefilter_pairs(
    const std::vector<GPElement>& catalog)
{
    std::vector<std::pair<uint32_t, uint32_t>> pairs;

    if (!config_.use_perigee_filter) {
        // No prefilter: return all pairs
        for (uint32_t i = 0; i < catalog.size(); i++) {
            for (uint32_t j = i + 1; j < catalog.size(); j++) {
                pairs.push_back({i, j});
            }
        }
        return pairs;
    }

    // Sort by perigee for sweep-line approach
    std::vector<uint32_t> indices(catalog.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
              [&catalog](uint32_t a, uint32_t b) {
                  return catalog[a].perigee_km < catalog[b].perigee_km;
              });

    // Sweep line: for each object, only check objects whose perigee
    // is within [perigee - margin, apogee + margin]
    const double margin = 50.0; // km
    for (size_t ii = 0; ii < indices.size(); ii++) {
        uint32_t i = indices[ii];
        double hi = catalog[i].apogee_km + margin;

        for (size_t jj = ii + 1; jj < indices.size(); jj++) {
            uint32_t j = indices[jj];
            if (catalog[j].perigee_km - margin > hi) break; // No more overlap

            if (altitude_overlap(catalog[i], catalog[j], margin)) {
                uint32_t a = std::min(i, j), b = std::max(i, j);
                pairs.push_back({a, b});
            }
        }
    }

    return pairs;
}

// ============================================================================
// Dynamic Windowing
// ============================================================================

double ConjunctionScreener::adaptive_step(double distance_km, double closing_rate_kms) const {
    if (!config_.use_dynamic_window) return config_.coarse_step_sec;

    // When objects are far apart and separating, use max step
    if (distance_km > config_.close_threshold_km * 2 || closing_rate_kms <= 0) {
        return config_.max_step_sec;
    }

    // When closing fast and nearby, reduce step
    // Step = max(min_step, distance / closing_rate / safety_factor)
    // This ensures we don't skip over the TCA
    if (closing_rate_kms > 0.01 && distance_km < config_.close_threshold_km) {
        double time_to_close = distance_km / closing_rate_kms; // seconds
        double step = time_to_close / 10.0; // Safety: 10 samples before potential TCA
        return std::max(config_.min_step_sec, std::min(step, config_.max_step_sec));
    }

    // Linear interpolation between min and max step based on distance
    double frac = std::min(1.0, distance_km / config_.close_threshold_km);
    return config_.min_step_sec + frac * (config_.max_step_sec - config_.min_step_sec);
}

// ============================================================================
// Worker Thread
// ============================================================================

ConjunctionScreener::ThreadWork ConjunctionScreener::process_time_steps(
    const std::vector<TLE>& tles,
    const std::vector<std::pair<uint32_t, uint32_t>>& valid_pairs,
    double start_jd, double end_jd, double step_sec,
    int thread_id, int total_threads)
{
    ThreadWork work;

    // Build position cache for all objects at each time step
    double step_days = step_sec / 86400.0;
    int total_steps = static_cast<int>((end_jd - start_jd) / step_days);

    // This thread handles steps: thread_id, thread_id + total_threads, ...
    for (int step = thread_id; step <= total_steps; step += total_threads) {
        double jd = start_jd + step * step_days;

        // Propagate all objects
        std::vector<KDPoint> points;
        points.reserve(tles.size());

        for (uint32_t i = 0; i < tles.size(); i++) {
            try {
                auto sv = propagate_sgp4(tles[i], jd);
                points.push_back({sv.x, sv.y, sv.z, i, static_cast<uint32_t>(step)});
                work.propagations++;
            } catch (...) {
                // SGP4 error — skip this object at this time
            }
        }

        if (config_.use_kdtree && points.size() > 50) {
            // Build KD-tree and find close pairs
            KDTree tree;
            tree.build(points);

            // Use 200x threshold for coarse KD-tree screening
            // LEO relative speeds of ~10 km/s mean 60s steps = 600 km travel
            auto candidates = tree.find_close_pairs(points, config_.threshold_km * 200.0);
            for (auto& c : candidates) {
                work.candidates.push_back(c);
            }
        } else {
            // Brute force for small catalogs — check all pairs of propagated points
            for (size_t ii = 0; ii < points.size(); ii++) {
                for (size_t jj = ii + 1; jj < points.size(); jj++) {
                    if (points[ii].obj_index == points[jj].obj_index) continue;

                    double dx = points[ii].x - points[jj].x;
                    double dy = points[ii].y - points[jj].y;
                    double dz = points[ii].z - points[jj].z;
                    double dist = std::sqrt(dx*dx + dy*dy + dz*dz);

                    // Use 10x threshold for coarse screening
                    // (objects at LEO speeds can be 50+ km apart at 60s intervals
                    // but still have sub-km TCA between steps)
                    double coarse_threshold = config_.threshold_km * 200.0;
                    if (dist <= coarse_threshold) {
                        uint32_t a = std::min(points[ii].obj_index, points[jj].obj_index);
                        uint32_t b = std::max(points[ii].obj_index, points[jj].obj_index);
                        work.candidates.push_back({a, b, static_cast<uint32_t>(step), dist});
                    }
                }
            }
        }
    }

    return work;
}

// ============================================================================
// Main Screening
// ============================================================================

std::vector<ConjunctionEvent> ConjunctionScreener::screen(
    const std::vector<GPElement>& catalog,
    ProgressCallback progress)
{
    return screen(catalog, catalog, progress);
}

std::vector<ConjunctionEvent> ConjunctionScreener::screen(
    const std::vector<GPElement>& primaries,
    const std::vector<GPElement>& secondaries,
    ProgressCallback progress)
{
    auto t_start = std::chrono::high_resolution_clock::now();
    stats_ = {};

    // Merge catalogs (dedup by NORAD ID)
    std::vector<GPElement> catalog;
    std::unordered_set<int> seen;
    for (const auto& gp : primaries) {
        if (seen.insert(gp.norad_cat_id).second) catalog.push_back(gp);
    }
    for (const auto& gp : secondaries) {
        if (seen.insert(gp.norad_cat_id).second) catalog.push_back(gp);
    }

    stats_.total_objects = catalog.size();
    if (progress) progress(0.05, "Prefiltering pairs by altitude...");

    // Step 1: Prefilter by altitude
    auto valid_pairs = prefilter_pairs(catalog);
    stats_.pairs_prefiltered = (catalog.size() * (catalog.size() - 1) / 2) - valid_pairs.size();
    stats_.pairs_screened = valid_pairs.size();

    if (progress) progress(0.1, "Converting GP to TLE...");

    // Step 2: Convert GP → TLE
    std::vector<TLE> tles;
    tles.reserve(catalog.size());
    for (const auto& gp : catalog) {
        try {
            tles.push_back(gp_to_tle(gp));
        } catch (...) {
            // Create a dummy TLE that will be skipped
            TLE dummy;
            dummy.norad_cat_id = gp.norad_cat_id;
            tles.push_back(dummy);
        }
    }

    double start_jd = config_.start_jd;
    double end_jd = start_jd + config_.duration_days;

    if (progress) progress(0.15, "Propagating and screening...");

    // Step 3: Parallel propagation + KD-tree screening
    int num_threads = std::max(1, config_.num_threads);
    std::vector<std::thread> threads;
    std::vector<ThreadWork> thread_results(num_threads);

    for (int t = 0; t < num_threads; t++) {
        threads.emplace_back([&, t]() {
            thread_results[t] = process_time_steps(
                tles, valid_pairs,
                start_jd, end_jd, config_.coarse_step_sec,
                t, num_threads);
        });
    }

    for (auto& t : threads) t.join();

    // Merge candidates from all threads
    std::vector<CandidatePair> all_candidates;
    for (const auto& tw : thread_results) {
        all_candidates.insert(all_candidates.end(),
                            tw.candidates.begin(), tw.candidates.end());
        stats_.propagations += tw.propagations;
    }

    if (progress) progress(0.7, "Deduplicating candidates...");

    // Deduplicate candidate pairs (same pair may appear at multiple time steps)
    // Keep the one with minimum distance
    std::map<uint64_t, CandidatePair> unique_pairs;
    for (const auto& c : all_candidates) {
        uint64_t key = (static_cast<uint64_t>(c.obj1_index) << 32) | c.obj2_index;
        auto it = unique_pairs.find(key);
        if (it == unique_pairs.end() || c.distance_km < it->second.distance_km) {
            unique_pairs[key] = c;
        }
    }

    stats_.kdtree_candidates = unique_pairs.size();

    if (progress) progress(0.75, "Refining TCA...");

    // Step 4: Fine TCA refinement for each candidate pair
    std::vector<ConjunctionEvent> events;
    std::mutex events_mutex;

    std::vector<std::pair<uint64_t, CandidatePair>> pair_list(
        unique_pairs.begin(), unique_pairs.end());

    // Parallel TCA refinement
    auto refine_range = [&](size_t from, size_t to) {
        std::vector<ConjunctionEvent> local_events;
        for (size_t i = from; i < to; i++) {
            const auto& [key, cand] = pair_list[i];
            try {
                auto event = assess_conjunction(
                    tles[cand.obj1_index], tles[cand.obj2_index],
                    start_jd, config_.duration_days,
                    config_.combined_radius_m / 2.0,
                    config_.combined_radius_m / 2.0);

                if (event.min_range_km <= config_.threshold_km) {
                    local_events.push_back(event);
                }
            } catch (...) {}
        }

        std::lock_guard<std::mutex> lock(events_mutex);
        events.insert(events.end(), local_events.begin(), local_events.end());
    };

    threads.clear();
    size_t chunk = (pair_list.size() + num_threads - 1) / num_threads;
    for (int t = 0; t < num_threads; t++) {
        size_t from = t * chunk;
        size_t to = std::min(from + chunk, pair_list.size());
        if (from < to) {
            threads.emplace_back(refine_range, from, to);
        }
    }
    for (auto& t : threads) t.join();

    stats_.tca_refined = pair_list.size();
    stats_.conjunctions_found = events.size();

    // Sort by max probability descending
    std::sort(events.begin(), events.end(),
              [](const ConjunctionEvent& a, const ConjunctionEvent& b) {
                  return a.max_probability > b.max_probability;
              });

    auto t_end = std::chrono::high_resolution_clock::now();
    stats_.elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    if (progress) progress(1.0, "Done");

    return events;
}

} // namespace conjunction
