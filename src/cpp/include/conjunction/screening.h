#ifndef CONJUNCTION_SCREENING_H
#define CONJUNCTION_SCREENING_H

/**
 * High-Performance Conjunction Screening Engine
 *
 * Architecture:
 *   1. Perigee/apogee prefilter (O(n) per pair — eliminate impossible pairs)
 *   2. Coarse propagation at ~60s steps with KD-tree spatial indexing
 *   3. Dynamic windowing: adapt step size based on closing rate
 *   4. Fine TCA search with golden section for candidates within threshold
 *   5. Full conjunction assessment for confirmed close approaches
 *
 * Threading (pthreads):
 *   - Time steps are divided among threads
 *   - Each thread propagates all objects for its assigned time steps
 *   - KD-tree built per time step (owned by thread)
 *   - Candidate pairs merged with mutex-free per-thread collection
 *   - TCA refinement parallelized per candidate pair
 *
 * For WASM (Emscripten):
 *   - -s USE_PTHREADS=1
 *   - -s PTHREAD_POOL_SIZE=navigator.hardwareConcurrency
 *   - -s SHARED_MEMORY=1
 *   - SharedArrayBuffer required in browser (COOP/COEP headers)
 */

#include "conjunction_assessment.h"
#include "kdtree.h"
#include "gp_json.h"
#include <vector>
#include <functional>
#include <atomic>

namespace conjunction {

/// Screening configuration
struct ScreeningConfig {
    double start_jd = 0.0;          // Start of screening window
    double duration_days = 7.0;     // Screening duration
    double threshold_km = 5.0;      // Miss distance threshold
    double coarse_step_sec = 60.0;  // Base coarse step size
    double fine_tol_sec = 0.001;    // TCA refinement tolerance
    double combined_radius_m = 10.0; // Combined hard-body radius
    int num_threads = 4;            // Thread pool size
    bool use_kdtree = true;         // Use KD-tree (vs brute force)
    bool use_dynamic_window = true; // Adaptive step size
    bool use_perigee_filter = true; // Prefilter by altitude overlap

    // Dynamic windowing params
    double min_step_sec = 5.0;      // Minimum step when objects closing fast
    double max_step_sec = 120.0;    // Maximum step when objects far apart
    double close_threshold_km = 100.0; // Distance to start reducing step
};

/// Screening progress callback
using ProgressCallback = std::function<void(double fraction, const std::string& status)>;

/// Screening statistics
struct ScreeningStats {
    uint64_t total_objects = 0;
    uint64_t pairs_screened = 0;
    uint64_t pairs_prefiltered = 0;  // Eliminated by perigee/apogee
    uint64_t kdtree_candidates = 0;  // Candidates from KD-tree
    uint64_t tca_refined = 0;        // Pairs with TCA refinement
    uint64_t conjunctions_found = 0; // Final confirmed conjunctions
    uint64_t propagations = 0;       // Total SGP4 propagations
    double elapsed_ms = 0.0;         // Wall clock time
};

/// Can objects at these altitudes ever meet?
/// Returns false if perigee/apogee ranges don't overlap
bool altitude_overlap(const GPElement& gp1, const GPElement& gp2,
                      double margin_km = 50.0);

/// High-performance conjunction screener
class ConjunctionScreener {
public:
    explicit ConjunctionScreener(const ScreeningConfig& config)
        : config_(config) {}

    /// Screen a catalog of GP elements for conjunctions
    /// Returns confirmed conjunction events sorted by max probability (desc)
    std::vector<ConjunctionEvent> screen(
        const std::vector<GPElement>& catalog,
        ProgressCallback progress = nullptr);

    /// Screen primary vs secondary catalog
    std::vector<ConjunctionEvent> screen(
        const std::vector<GPElement>& primaries,
        const std::vector<GPElement>& secondaries,
        ProgressCallback progress = nullptr);

    /// Get statistics from the last screening run
    const ScreeningStats& stats() const { return stats_; }

private:
    ScreeningConfig config_;
    ScreeningStats stats_;
    std::atomic<uint64_t> propagation_count_{0};

    /// Prefilter pairs by perigee/apogee overlap
    std::vector<std::pair<uint32_t, uint32_t>> prefilter_pairs(
        const std::vector<GPElement>& catalog);

    /// Worker thread: propagate objects for assigned time steps,
    /// build KD-tree, find candidates
    struct ThreadWork {
        std::vector<CandidatePair> candidates;
        uint64_t propagations = 0;
    };

    ThreadWork process_time_steps(
        const std::vector<TLE>& tles,
        const std::vector<std::pair<uint32_t, uint32_t>>& valid_pairs,
        double start_jd, double end_jd, double step_sec,
        int thread_id, int total_threads);

    /// Dynamic step size based on closing rate
    double adaptive_step(double distance_km, double closing_rate_kms) const;
};

} // namespace conjunction

#endif // CONJUNCTION_SCREENING_H
