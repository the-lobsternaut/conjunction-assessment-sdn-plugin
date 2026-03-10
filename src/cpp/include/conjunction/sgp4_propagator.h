#ifndef CONJUNCTION_SGP4_PROPAGATOR_H
#define CONJUNCTION_SGP4_PROPAGATOR_H

/**
 * SGP4 Propagator — Lightweight wrapper around Tudat's TLE/SGP4 implementation.
 *
 * Source: DigitalArsenal/tudat (fork of tudat-team/tudat)
 * Files:
 *   - include/tudat/astro/ephemerides/tleEphemeris.h
 *   - src/astro/ephemerides/tleEphemeris.cpp
 *   - include/tudat/astro/basic_astro/tleElementsConversions.h
 *
 * For WASM build, we extract only the SGP4 propagation code from Tudat
 * to avoid pulling in the full Tudat dependency tree.
 */

#include <string>
#include <vector>
#include <cstdint>

namespace conjunction {

/// State vector: position (km) + velocity (km/s) in J2000/TEME
struct StateVector {
    double epoch_jd;   // Julian Date
    double x, y, z;    // Position (km)
    double vx, vy, vz; // Velocity (km/s)
};

/// Two-Line Element set
struct TLE {
    std::string name;
    std::string line1;
    std::string line2;

    // Parsed fields
    int norad_cat_id = 0;
    double epoch_jd = 0.0;      // Julian Date of epoch
    double bstar = 0.0;
    double inclination = 0.0;   // degrees
    double raan = 0.0;          // degrees
    double eccentricity = 0.0;
    double arg_perigee = 0.0;   // degrees
    double mean_anomaly = 0.0;  // degrees
    double mean_motion = 0.0;   // rev/day
};

// Forward declaration
struct GPElement;

/// Parse a TLE from 3 lines (name + line1 + line2)
TLE parse_tle(const std::string& name, const std::string& line1, const std::string& line2);

/// Parse all TLEs from a multi-line string
std::vector<TLE> parse_tle_file(const std::string& data);

/// Propagate a TLE to a given Julian Date using SGP4
/// Returns state in TEME frame (km, km/s)
StateVector propagate_sgp4(const TLE& tle, double target_jd);

/// Propagate directly from GP/OMM elements (no TLE text intermediary)
/// Supports any NORAD catalog ID (integer, no 5-digit limit)
StateVector propagate_sgp4_gp(const GPElement& gp, double target_jd);

/// Convert Julian Date to ISO 8601 string
std::string jd_to_iso(double jd);

/// Convert ISO 8601 string to Julian Date
double iso_to_jd(const std::string& iso);

/// Get Julian Date for a year + fractional day of year
double epoch_to_jd(int year, double day_of_year);

} // namespace conjunction

#endif // CONJUNCTION_SGP4_PROPAGATOR_H
