/**
 * EphemerisSource — OEM interpolation implementation
 */

#include "conjunction/ephemeris_source.h"
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace conjunction {

// ── OEM Hermite Interpolation ──

StateVector OEMEphemerisSource::state_at(double jd) const {
    if (points_.empty()) {
        throw std::runtime_error("Empty ephemeris");
    }

    // Clamp to valid range (extrapolate to nearest point)
    if (jd <= points_.front().jd) {
        const auto& p = points_.front();
        return {p.jd, p.x, p.y, p.z, p.vx, p.vy, p.vz};
    }
    if (jd >= points_.back().jd) {
        const auto& p = points_.back();
        return {p.jd, p.x, p.y, p.z, p.vx, p.vy, p.vz};
    }

    // Binary search for bracketing interval
    size_t lo = 0, hi = points_.size() - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (points_[mid].jd <= jd) lo = mid;
        else hi = mid;
    }

    // Hermite interpolation between points_[lo] and points_[hi]
    const auto& p0 = points_[lo];
    const auto& p1 = points_[hi];

    double dt = p1.jd - p0.jd;
    if (dt < 1e-15) {
        return {jd, p0.x, p0.y, p0.z, p0.vx, p0.vy, p0.vz};
    }

    double t = (jd - p0.jd) / dt;  // normalized [0,1]
    double t2 = t * t;
    double t3 = t2 * t;

    // Hermite basis functions
    double h00 = 2*t3 - 3*t2 + 1;
    double h10 = t3 - 2*t2 + t;
    double h01 = -2*t3 + 3*t2;
    double h11 = t3 - t2;

    // dt in seconds for velocity scaling
    double dt_sec = dt * 86400.0;

    StateVector sv;
    sv.epoch_jd = jd;

    // Position interpolation
    sv.x = h00 * p0.x + h10 * p0.vx * dt_sec + h01 * p1.x + h11 * p1.vx * dt_sec;
    sv.y = h00 * p0.y + h10 * p0.vy * dt_sec + h01 * p1.y + h11 * p1.vy * dt_sec;
    sv.z = h00 * p0.z + h10 * p0.vz * dt_sec + h01 * p1.z + h11 * p1.vz * dt_sec;

    // Velocity from Hermite derivative
    double dh00 = 6*t2 - 6*t;
    double dh10 = 3*t2 - 4*t + 1;
    double dh01 = -6*t2 + 6*t;
    double dh11 = 3*t2 - 2*t;

    sv.vx = (dh00 * p0.x + dh10 * p0.vx * dt_sec + dh01 * p1.x + dh11 * p1.vx * dt_sec) / dt_sec;
    sv.vy = (dh00 * p0.y + dh10 * p0.vy * dt_sec + dh01 * p1.y + dh11 * p1.vy * dt_sec) / dt_sec;
    sv.vz = (dh00 * p0.z + dh10 * p0.vz * dt_sec + dh01 * p1.z + dh11 * p1.vz * dt_sec) / dt_sec;

    return sv;
}

} // namespace conjunction
