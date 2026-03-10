#ifndef CONJUNCTION_KDTREE_H
#define CONJUNCTION_KDTREE_H

/**
 * KD-Tree for Spatial Conjunction Screening
 *
 * Uses a 3D KD-tree to find potential conjunctions in O(n log n) average
 * case instead of O(n²) brute force. The tree is built from propagated
 * positions at discrete time steps and queried with a range threshold.
 *
 * Key optimization: Objects are inserted with their positions at each
 * time step. The tree is rebuilt at each time step (or a sliding window
 * approach is used). Pairs found within the screening threshold are
 * then refined with fine-grained TCA search.
 *
 * Thread Safety:
 *   - Tree build: thread-safe (each thread builds own tree, or use
 *     shared read-only tree with write barrier between steps)
 *   - Range query: read-only, fully thread-safe
 *   - Pair refinement: each pair processed independently (embarrassingly parallel)
 *
 * Compile with -pthread for pthreads support.
 * For WASM: -s USE_PTHREADS=1 -s PTHREAD_POOL_SIZE=N
 */

#include <vector>
#include <cstdint>
#include <cstddef>
#include <algorithm>
#include <cmath>

namespace conjunction {

/// 3D point with object index metadata
struct KDPoint {
    double x, y, z;      // Position (km)
    uint32_t obj_index;   // Index into the object array
    uint32_t time_step;   // Time step index
};

/// KD-Tree node (implicit array layout for cache efficiency)
struct KDNode {
    double split_val;     // Split value on split_axis
    uint32_t left;        // Left child index (0 = no child)
    uint32_t right;       // Right child index (0 = no child)
    uint32_t point_idx;   // Index into points array
    uint8_t split_axis;   // 0=x, 1=y, 2=z
};

/// Candidate pair from KD-tree screening
struct CandidatePair {
    uint32_t obj1_index;
    uint32_t obj2_index;
    uint32_t time_step;
    double distance_km;   // Approximate distance at this time step
};

/// KD-Tree for spatial range queries
class KDTree {
public:
    KDTree() = default;

    /// Build tree from a set of points
    /// Points are sorted in-place during construction
    void build(std::vector<KDPoint>& points);

    /// Find all points within `radius_km` of `query`
    /// Results are appended to `results`
    void range_query(const KDPoint& query, double radius_km,
                     std::vector<uint32_t>& results) const;

    /// Find all pairs of points within `radius_km` of each other
    /// This is the main screening function — returns candidate pairs
    /// sorted by distance
    std::vector<CandidatePair> find_close_pairs(
        const std::vector<KDPoint>& points,
        double radius_km) const;

    /// Clear the tree
    void clear() { nodes_.clear(); }

    /// Number of nodes
    size_t size() const { return nodes_.size(); }

private:
    std::vector<KDNode> nodes_;
    std::vector<KDPoint> sorted_points_;

    uint32_t build_recursive(std::vector<KDPoint>& points,
                            size_t start, size_t end, int depth);

    void range_query_recursive(uint32_t node_idx, const KDPoint& query,
                              double radius_sq,
                              std::vector<uint32_t>& results) const;

    static double point_distance_sq(const KDPoint& a, const KDPoint& b) {
        double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
        return dx*dx + dy*dy + dz*dz;
    }

    static double point_axis(const KDPoint& p, int axis) {
        switch (axis) {
            case 0: return p.x;
            case 1: return p.y;
            case 2: return p.z;
            default: return 0;
        }
    }
};

} // namespace conjunction

#endif // CONJUNCTION_KDTREE_H
