#ifndef CONJUNCTION_GP_JSON_H
#define CONJUNCTION_GP_JSON_H

/**
 * GP JSON Parser — Parse CelesTrak GP data in JSON/CSV format
 *
 * CelesTrak now provides GP data as structured JSON:
 *   https://celestrak.org/NORAD/elements/gp.php?CATNR=25544&FORMAT=json
 *
 * Fields:
 *   OBJECT_NAME, OBJECT_ID, EPOCH, MEAN_MOTION, ECCENTRICITY,
 *   INCLINATION, RA_OF_ASC_NODE, ARG_OF_PERICENTER, MEAN_ANOMALY,
 *   EPHEMERIS_TYPE, CLASSIFICATION_TYPE, NORAD_CAT_ID, ELEMENT_SET_NO,
 *   REV_AT_EPOCH, BSTAR, MEAN_MOTION_DOT, MEAN_MOTION_DDOT
 *
 * This eliminates the need for TLE text line parsing and CelesTrak's
 * non-standard 70-char line normalization.
 */

#include "sgp4_propagator.h"
#include <string>
#include <vector>

namespace conjunction {

/// GP element set (structured, mirrors CelesTrak JSON fields)
struct GPElement {
    std::string object_name;
    std::string object_id;       // International designator
    std::string epoch_iso;       // ISO 8601 epoch
    double mean_motion = 0.0;    // rev/day
    double eccentricity = 0.0;
    double inclination = 0.0;    // degrees
    double ra_of_asc_node = 0.0; // degrees
    double arg_of_pericenter = 0.0; // degrees
    double mean_anomaly = 0.0;   // degrees
    int ephemeris_type = 0;
    char classification_type = 'U';
    int norad_cat_id = 0;
    int element_set_no = 0;
    int rev_at_epoch = 0;
    double bstar = 0.0;
    double mean_motion_dot = 0.0;
    double mean_motion_ddot = 0.0;

    // Derived fields
    double epoch_jd = 0.0;
    double semi_major_axis_km = 0.0;  // Computed from mean motion
    double perigee_km = 0.0;          // For KD-tree prefilter
    double apogee_km = 0.0;           // For KD-tree prefilter
};

/// Parse GP JSON array (CelesTrak format)
/// Input: "[{...}, {...}, ...]"
std::vector<GPElement> parse_gp_json(const std::string& json);

/// Parse GP CSV (CelesTrak format)
/// Input: "OBJECT_NAME,OBJECT_ID,...\nISS,1998-067A,..."
std::vector<GPElement> parse_gp_csv(const std::string& csv);

/// Convert GPElement to TLE struct (for SGP4 propagation via dnwrnr/sgp4)
/// Reconstructs TLE line1 + line2 text from structured GP fields
TLE gp_to_tle(const GPElement& gp);

/// Compute derived orbital parameters from GP elements
/// (semi-major axis, perigee, apogee altitudes)
void compute_derived(GPElement& gp);

} // namespace conjunction

#endif // CONJUNCTION_GP_JSON_H
