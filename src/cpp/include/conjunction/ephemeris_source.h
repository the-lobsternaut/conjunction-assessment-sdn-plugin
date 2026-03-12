#ifndef CONJUNCTION_EPHEMERIS_SOURCE_H
#define CONJUNCTION_EPHEMERIS_SOURCE_H

/**
 * EphemerisSource — Propagator-agnostic state provider
 *
 * Any propagator (SGP4, HPOP, tabulated OEM, external API) implements
 * this interface. The conjunction assessment engine only sees states.
 *
 * Implementations:
 *   SGP4EphemerisSource     — from TLE/GP elements
 *   OEMEphemerisSource      — from OEM ephemeris (interpolated)
 *   StateEphemerisSource    — from pre-computed state vector (two-body)
 *   CallbackEphemerisSource — from user-supplied function pointer
 */

#include "conjunction/sgp4_propagator.h"  // for StateVector
#include "conjunction/gp_json.h"          // for GPElement
#include <memory>
#include <string>
#include <vector>
#include <functional>

namespace conjunction {

// ── Abstract Interface ──

class EphemerisSource {
public:
    virtual ~EphemerisSource() = default;

    /// Get state at Julian Date. Returns position (km) + velocity (km/s).
    virtual StateVector state_at(double jd) const = 0;

    /// Epoch JD of the underlying data (for DSE computation)
    virtual double epoch_jd() const = 0;

    /// Object identifier (NORAD ID, name, etc.)
    virtual std::string object_id() const { return ""; }
    virtual std::string object_name() const { return ""; }
    virtual int norad_id() const { return 0; }

    /// Valid time range [start_jd, end_jd]. 0 = unlimited.
    virtual double valid_start_jd() const { return 0; }
    virtual double valid_end_jd() const { return 0; }
};

// ── SGP4 from TLE ──

class SGP4EphemerisSource : public EphemerisSource {
public:
    explicit SGP4EphemerisSource(const TLE& tle) : tle_(tle) {}

    StateVector state_at(double jd) const override {
        return propagate_sgp4(tle_, jd);
    }

    double epoch_jd() const override { return tle_.epoch_jd; }
    std::string object_id() const override { return std::to_string(tle_.norad_cat_id); }
    std::string object_name() const override { return tle_.name; }
    int norad_id() const override { return tle_.norad_cat_id; }

private:
    TLE tle_;
};

// ── SGP4 from GP/OMM elements ──

class GPEphemerisSource : public EphemerisSource {
public:
    explicit GPEphemerisSource(const GPElement& gp) : gp_(gp) {}

    StateVector state_at(double jd) const override {
        return propagate_sgp4_gp(gp_, jd);
    }

    double epoch_jd() const override { return gp_.epoch_jd; }
    std::string object_id() const override { return std::to_string(gp_.norad_cat_id); }
    std::string object_name() const override { return gp_.object_name; }
    int norad_id() const override { return gp_.norad_cat_id; }

private:
    GPElement gp_;
};

// ── Tabulated OEM (Hermite interpolation) ──

struct EphemerisPoint {
    double jd;
    double x, y, z;
    double vx, vy, vz;
};

class OEMEphemerisSource : public EphemerisSource {
public:
    OEMEphemerisSource(std::vector<EphemerisPoint> points,
                        std::string name = "", int norad = 0)
        : points_(std::move(points)), name_(std::move(name)), norad_(norad) {
        if (!points_.empty()) {
            epoch_jd_ = points_.front().jd;
        }
    }

    StateVector state_at(double jd) const override;

    double epoch_jd() const override { return epoch_jd_; }
    std::string object_name() const override { return name_; }
    int norad_id() const override { return norad_; }
    double valid_start_jd() const override {
        return points_.empty() ? 0 : points_.front().jd;
    }
    double valid_end_jd() const override {
        return points_.empty() ? 0 : points_.back().jd;
    }

private:
    std::vector<EphemerisPoint> points_;
    double epoch_jd_ = 0;
    std::string name_;
    int norad_ = 0;
};

// ── Callback-based (user-supplied function) ──

using EphemerisCallback = std::function<StateVector(double jd)>;

class CallbackEphemerisSource : public EphemerisSource {
public:
    CallbackEphemerisSource(EphemerisCallback cb, double epoch = 0,
                             std::string name = "", int norad = 0)
        : cb_(std::move(cb)), epoch_(epoch), name_(std::move(name)), norad_(norad) {}

    StateVector state_at(double jd) const override { return cb_(jd); }
    double epoch_jd() const override { return epoch_; }
    std::string object_name() const override { return name_; }
    int norad_id() const override { return norad_; }

private:
    EphemerisCallback cb_;
    double epoch_;
    std::string name_;
    int norad_;
};

} // namespace conjunction

#endif // CONJUNCTION_EPHEMERIS_SOURCE_H
