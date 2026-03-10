/**
 * KD-Tree for Spatial Conjunction Screening
 *
 * 3D KD-tree optimized for satellite position queries.
 * Uses implicit array layout for cache efficiency.
 * All operations are read-only after build — thread-safe for queries.
 */

#include "conjunction/kdtree.h"
#include <algorithm>
#include <numeric>

namespace conjunction {

void KDTree::build(std::vector<KDPoint>& points) {
    nodes_.clear();

    if (points.empty()) return;

    // Reserve space (binary tree: at most 2N nodes)
    nodes_.reserve(points.size() * 2);

    build_recursive(points, 0, points.size(), 0);

    // Store the sorted points (modified in place by nth_element)
    sorted_points_ = points;
}

uint32_t KDTree::build_recursive(std::vector<KDPoint>& points,
                                 size_t start, size_t end, int depth) {
    if (start >= end) return 0;

    int axis = depth % 3;
    size_t mid = (start + end) / 2;

    // Partial sort to find median on current axis
    std::nth_element(points.begin() + start, points.begin() + mid,
                     points.begin() + end,
                     [axis](const KDPoint& a, const KDPoint& b) {
                         return point_axis(a, axis) < point_axis(b, axis);
                     });

    // Create node
    uint32_t node_idx = static_cast<uint32_t>(nodes_.size());
    nodes_.push_back({});
    KDNode& node = nodes_[node_idx];

    node.split_axis = axis;
    node.split_val = point_axis(points[mid], axis);
    node.point_idx = static_cast<uint32_t>(mid);

    // Recurse
    // Need to handle carefully since nodes_ may reallocate
    uint32_t left_idx = 0, right_idx = 0;

    if (mid > start) {
        left_idx = build_recursive(points, start, mid, depth + 1);
    }
    if (mid + 1 < end) {
        right_idx = build_recursive(points, mid + 1, end, depth + 1);
    }

    // Set children (after recursion, since nodes_ may have reallocated)
    nodes_[node_idx].left = left_idx;
    nodes_[node_idx].right = right_idx;

    return node_idx;
}

void KDTree::range_query(const KDPoint& query, double radius_km,
                         std::vector<uint32_t>& results) const {
    if (nodes_.empty()) return;
    double radius_sq = radius_km * radius_km;
    range_query_recursive(0, query, radius_sq, results);
}

void KDTree::range_query_recursive(uint32_t node_idx, const KDPoint& query,
                                   double radius_sq,
                                   std::vector<uint32_t>& results) const {
    if (node_idx >= nodes_.size()) return;

    const KDNode& node = nodes_[node_idx];
    const KDPoint& point = sorted_points_[node.point_idx];

    // Check this point
    double dist_sq = point_distance_sq(query, point);
    if (dist_sq <= radius_sq) {
        results.push_back(node.point_idx);
    }

    // Which side of the split are we on?
    double query_val = point_axis(query, node.split_axis);
    double diff = query_val - node.split_val;

    // Visit closer child first
    uint32_t near_child = (diff < 0) ? node.left : node.right;
    uint32_t far_child = (diff < 0) ? node.right : node.left;

    if (near_child != 0) {
        range_query_recursive(near_child, query, radius_sq, results);
    }

    // Visit far child only if splitting plane is within radius
    if (far_child != 0 && diff * diff <= radius_sq) {
        range_query_recursive(far_child, query, radius_sq, results);
    }
}

std::vector<CandidatePair> KDTree::find_close_pairs(
    const std::vector<KDPoint>& /*points*/,
    double radius_km) const
{
    std::vector<CandidatePair> pairs;
    if (sorted_points_.empty()) return pairs;

    // For each point in the tree, query for neighbors within radius
    // Skip self-matches and duplicate pairs (only keep obj_i < obj_j)
    std::vector<uint32_t> neighbors;

    for (size_t i = 0; i < sorted_points_.size(); i++) {
        const KDPoint& pt = sorted_points_[i];
        neighbors.clear();
        range_query(pt, radius_km, neighbors);

        for (uint32_t nb_idx : neighbors) {
            const KDPoint& neighbor = sorted_points_[nb_idx];

            // Skip self
            if (neighbor.obj_index == pt.obj_index) continue;

            // Only keep pair once (lower index first)
            if (pt.obj_index >= neighbor.obj_index) continue;

            double dist = std::sqrt(point_distance_sq(pt, neighbor));

            pairs.push_back({
                pt.obj_index,
                neighbor.obj_index,
                pt.time_step,
                dist
            });
        }
    }

    // Sort by distance
    std::sort(pairs.begin(), pairs.end(),
              [](const CandidatePair& a, const CandidatePair& b) {
                  return a.distance_km < b.distance_km;
              });

    return pairs;
}

} // namespace conjunction
