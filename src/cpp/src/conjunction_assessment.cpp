/**
 * Conjunction Assessment Engine Implementation
 *
 * Uses dnwrnr/sgp4 (Apache 2.0) for TLE parsing + SGP4 propagation.
 * Implements Alfano maximum probability method (AAS 03-548).
 * Validates against CelesTrak SOCRATES Plus.
 */

#include "conjunction/conjunction_assessment.h"

// dnwrnr/sgp4 headers
#include "Tle.h"
#include "SGP4.h"
#include "DateTime.h"
#include "TimeSpan.h"
#include "Eci.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <iomanip>

namespace conjunction {

// ============================================================================
// Time Utilities
// ============================================================================

double epoch_to_jd(int year, double day_of_year) {
    // Julian Date from year + day of year
    int y = year;
    if (y < 57) y += 2000; else if (y < 100) y += 1900;

    int a = (14 - 1) / 12;
    int yy = y + 4800 - a;
    int mm = 1 + 12 * a - 3;
    double jd = 1 + (153 * mm + 2) / 5 + 365 * yy + yy / 4 - yy / 100 + yy / 400 - 32045;
    return jd - 0.5 + (day_of_year - 1.0);
}

std::string jd_to_iso(double jd) {
    // JD → calendar date
    double z = std::floor(jd + 0.5);
    double f = (jd + 0.5) - z;
    double a;
    if (z < 2299161) { a = z; }
    else {
        double alpha = std::floor((z - 1867216.25) / 36524.25);
        a = z + 1 + alpha - std::floor(alpha / 4);
    }
    double b = a + 1524;
    double c = std::floor((b - 122.1) / 365.25);
    double d = std::floor(365.25 * c);
    double e = std::floor((b - d) / 30.6001);

    double day = b - d - std::floor(30.6001 * e) + f;
    int month = (e < 14) ? (int)e - 1 : (int)e - 13;
    int year = (month > 2) ? (int)c - 4716 : (int)c - 4715;

    int day_int = (int)day;
    double frac = day - day_int;
    int hour = (int)(frac * 24);
    int minute = (int)((frac * 24 - hour) * 60);
    double second = ((frac * 24 - hour) * 60 - minute) * 60;

    std::stringstream ss;
    ss << std::setfill('0')
       << std::setw(4) << year << "-"
       << std::setw(2) << month << "-"
       << std::setw(2) << day_int << "T"
       << std::setw(2) << hour << ":"
       << std::setw(2) << minute << ":"
       << std::fixed << std::setprecision(3) << std::setw(6) << std::setfill('0') << second
       << "Z";
    return ss.str();
}

double iso_to_jd(const std::string& iso) {
    int year, month, day, hour, minute;
    double second;
    sscanf(iso.c_str(), "%d-%d-%dT%d:%d:%lf", &year, &month, &day, &hour, &minute, &second);

    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;
    double jd = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
    return jd - 0.5 + (hour + minute / 60.0 + second / 3600.0) / 24.0;
}

// ============================================================================
// TLE Parsing (wraps dnwrnr/sgp4)
// ============================================================================

/// Normalize a TLE line to valid format (exactly 69 chars)
static std::string normalize_tle_line(const std::string& line) {
    std::string l = line;
    // Trim trailing whitespace/CR
    while (!l.empty() && (l.back() == ' ' || l.back() == '\r' || l.back() == '\n'))
        l.pop_back();
    // Trim leading whitespace
    size_t start = 0;
    while (start < l.size() && (l[start] == ' ' || l[start] == '\t'))
        start++;
    if (start > 0) l = l.substr(start);
    // Handle 70-char lines from CelesTrak
    // Some TLEs have 8-digit eccentricity (non-standard) making line 70 chars.
    // The standard eccentricity field is 7 digits at columns 26-32.
    // If line is 70 chars, check for 8-digit eccentricity and truncate to 7.
    if (l.length() == 70 && l[0] == '2') {
        // Check if cols 26-33 are all digits (8-digit eccentricity)
        bool all_digits = true;
        for (int k = 26; k < 34 && k < (int)l.length(); k++) {
            if (!std::isdigit(l[k])) { all_digits = false; break; }
        }
        if (all_digits) {
            // Remove the last digit of the eccentricity field (col 33)
            l.erase(33, 1);
        } else {
            // Try removing a trailing space in the rev number area
            size_t pos = l.rfind("  ");
            if (pos >= 60) l.erase(pos, 1);
        }
    }
    if (l.length() == 70 && l[0] == '1') {
        size_t pos = l.rfind("  ");
        if (pos >= 60) l.erase(pos, 1);
    }
    // Pad to 69 chars if too short
    if (l.length() < 69) l.resize(69, ' ');
    // Truncate if still too long
    if (l.length() > 69) l = l.substr(0, 69);
    return l;
}

TLE parse_tle(const std::string& name, const std::string& line1, const std::string& line2) {
    TLE tle;
    tle.name = name;
    tle.line1 = normalize_tle_line(line1);
    tle.line2 = normalize_tle_line(line2);

    // Use dnwrnr Tle parser for the heavy lifting
    try {
        libsgp4::Tle sgp4_tle(name, tle.line1, tle.line2);
        tle.norad_cat_id = sgp4_tle.NoradNumber();
        tle.inclination = sgp4_tle.Inclination(true);   // degrees
        tle.raan = sgp4_tle.RightAscendingNode(true);
        tle.eccentricity = sgp4_tle.Eccentricity();
        tle.arg_perigee = sgp4_tle.ArgumentPerigee(true);
        tle.mean_anomaly = sgp4_tle.MeanAnomaly(true);
        tle.mean_motion = sgp4_tle.MeanMotion();         // rev/day
        tle.bstar = sgp4_tle.BStar();

        // Compute epoch JD
        auto epoch = sgp4_tle.Epoch();
        tle.epoch_jd = epoch.ToJulian();
    } catch (const std::exception& e) {
        throw std::runtime_error("TLE parse error: " + std::string(e.what()));
    }

    return tle;
}

std::vector<TLE> parse_tle_file(const std::string& data) {
    std::vector<TLE> tles;
    std::istringstream stream(data);
    std::string line;
    std::vector<std::string> lines;

    while (std::getline(stream, line)) {
        // Trim trailing whitespace/CR
        while (!line.empty() && (line.back() == '\r' || line.back() == ' ' || line.back() == '\n'))
            line.pop_back();
        // Trim leading whitespace (some GP data has leading spaces on line2)
        size_t start = 0;
        while (start < line.size() && (line[start] == ' ' || line[start] == '\t'))
            start++;
        if (start > 0) line = line.substr(start);
        if (!line.empty())
            lines.push_back(line);
    }

    // Parse as 3-line sets (name + line1 + line2)
    for (size_t i = 0; i + 2 < lines.size(); ) {
        if (lines[i + 1].size() > 1 && lines[i + 1][0] == '1' &&
            lines[i + 2].size() > 1 && lines[i + 2][0] == '2') {
            try {
                tles.push_back(parse_tle(lines[i], lines[i + 1], lines[i + 2]));
            } catch (...) {
                // Skip malformed TLEs
            }
            i += 3;
        } else if (lines[i].size() > 1 && lines[i][0] == '1' &&
                   lines[i + 1].size() > 1 && lines[i + 1][0] == '2') {
            // 2-line format (no name)
            try {
                tles.push_back(parse_tle("", lines[i], lines[i + 1]));
            } catch (...) {}
            i += 2;
        } else {
            i++;
        }
    }

    return tles;
}

// ============================================================================
// SGP4 Propagation
// ============================================================================

StateVector propagate_sgp4(const TLE& tle, double target_jd) {
    libsgp4::Tle sgp4_tle(tle.name, tle.line1, tle.line2);
    libsgp4::SGP4 sgp4(sgp4_tle);

    // Compute time since TLE epoch
    double delta_days = target_jd - tle.epoch_jd;
    int whole_days = static_cast<int>(delta_days);
    double frac_day = delta_days - whole_days;
    int hours = static_cast<int>(frac_day * 24.0);
    double frac_hour = frac_day * 24.0 - hours;
    int minutes = static_cast<int>(frac_hour * 60.0);
    double frac_min = frac_hour * 60.0 - minutes;
    int seconds = static_cast<int>(frac_min * 60.0);
    int microseconds = static_cast<int>((frac_min * 60.0 - seconds) * 1000000.0);

    libsgp4::DateTime epoch = sgp4_tle.Epoch();
    libsgp4::TimeSpan offset(whole_days, hours, minutes, seconds, microseconds);
    libsgp4::DateTime target_dt = epoch + offset;

    // Propagate
    libsgp4::Eci eci = sgp4.FindPosition(target_dt);

    StateVector sv;
    sv.epoch_jd = target_jd;
    sv.x = eci.Position().x;
    sv.y = eci.Position().y;
    sv.z = eci.Position().z;
    sv.vx = eci.Velocity().x;
    sv.vy = eci.Velocity().y;
    sv.vz = eci.Velocity().z;
    return sv;
}

// ============================================================================
// Frame Transformations
// ============================================================================

static void cross(const double a[3], const double b[3], double c[3]) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

static double dot(const double a[3], const double b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static double norm(const double a[3]) {
    return std::sqrt(dot(a, a));
}

static void normalize(double a[3]) {
    double n = norm(a);
    if (n > 1e-15) { a[0] /= n; a[1] /= n; a[2] /= n; }
}

void inertial_to_rtn(const StateVector& ref, const StateVector& target,
                     double& r, double& t, double& n_out,
                     double& vr, double& vt, double& vn) {
    // R = radial (along position vector)
    // T = transverse (along velocity, ⊥ R)
    // N = normal (R × T)
    double R_hat[3] = { ref.x, ref.y, ref.z };
    normalize(R_hat);

    double h[3]; // angular momentum
    double rv[3] = { ref.x, ref.y, ref.z };
    double vv[3] = { ref.vx, ref.vy, ref.vz };
    cross(rv, vv, h);

    double N_hat[3] = { h[0], h[1], h[2] };
    normalize(N_hat);

    double T_hat[3];
    cross(N_hat, R_hat, T_hat);
    normalize(T_hat);

    // Relative position
    double dr[3] = {
        target.x - ref.x,
        target.y - ref.y,
        target.z - ref.z
    };

    r = dot(dr, R_hat);
    t = dot(dr, T_hat);
    n_out = dot(dr, N_hat);

    // Relative velocity
    double dv[3] = {
        target.vx - ref.vx,
        target.vy - ref.vy,
        target.vz - ref.vz
    };

    vr = dot(dv, R_hat);
    vt = dot(dv, T_hat);
    vn = dot(dv, N_hat);
}

// ============================================================================
// TCA Finding
// ============================================================================

static double distance_at_jd(const TLE& tle1, const TLE& tle2, double jd) {
    auto s1 = propagate_sgp4(tle1, jd);
    auto s2 = propagate_sgp4(tle2, jd);
    double dx = s1.x - s2.x;
    double dy = s1.y - s2.y;
    double dz = s1.z - s2.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/// Find all local distance minima (potential close approaches)
static std::vector<std::pair<double, double>> find_all_minima(
    const TLE& tle1, const TLE& tle2,
    double start_jd, double end_jd, double step_sec)
{
    double step = step_sec / 86400.0;
    std::vector<std::pair<double, double>> minima; // (jd, distance)

    double prev_d = 1e18, curr_d = 1e18, prev_jd = start_jd;
    bool prev_decreasing = false;

    for (double jd = start_jd; jd <= end_jd; jd += step) {
        try {
            double d = distance_at_jd(tle1, tle2, jd);
            bool decreasing = (d < curr_d);

            // Local minimum: was decreasing, now increasing
            if (prev_decreasing && !decreasing && curr_d < 50.0) {
                minima.push_back({prev_jd, curr_d});
            }

            prev_decreasing = decreasing;
            prev_d = curr_d;
            curr_d = d;
            prev_jd = jd;
        } catch (...) {}
    }

    return minima;
}

/// Refine a local minimum using golden section search
static double refine_minimum(const TLE& tle1, const TLE& tle2,
                            double center_jd, double window_days,
                            double tol_sec)
{
    double a = center_jd - window_days;
    double b = center_jd + window_days;
    double tol = tol_sec / 86400.0;
    const double phi = (std::sqrt(5.0) - 1.0) / 2.0;

    while ((b - a) > tol) {
        double c = b - phi * (b - a);
        double d_val = a + phi * (b - a);

        try {
            double fc = distance_at_jd(tle1, tle2, c);
            double fd = distance_at_jd(tle1, tle2, d_val);

            if (fc < fd) b = d_val;
            else a = c;
        } catch (...) { break; }
    }

    return (a + b) / 2.0;
}

double find_tca(const TLE& tle1, const TLE& tle2,
                double start_jd, double duration_days,
                double coarse_step_sec, double fine_tol_sec) {
    double end_jd = start_jd + duration_days;
    double step = coarse_step_sec;

    // Step 1: Find all local minima with 10-second steps
    // (LEO orbital period ~90 min, so 10s gives ~540 samples per orbit)
    auto minima = find_all_minima(tle1, tle2, start_jd, end_jd, 10.0);

    if (minima.empty()) {
        // Fall back to global minimum with coarser step
        minima = find_all_minima(tle1, tle2, start_jd, end_jd, step);
    }

    if (minima.empty()) {
        return start_jd; // No close approach found
    }

    // Step 2: Find the closest approach among all minima
    double best_jd = minima[0].first;
    double best_dist = minima[0].second;

    for (const auto& m : minima) {
        if (m.second < best_dist) {
            best_dist = m.second;
            best_jd = m.first;
        }
    }

    // Step 3: Refine with golden section search
    double step_days = 10.0 / 86400.0;
    return refine_minimum(tle1, tle2, best_jd, step_days * 2, fine_tol_sec);
}

// ============================================================================
// Alfano Maximum Probability (AAS 03-548)
// ============================================================================

ProbResult alfano_max_probability(double miss_distance_km, double combined_radius_km) {
    ProbResult result;

    double d = miss_distance_km * 1000.0;  // Convert to meters
    double Rc = combined_radius_km * 1000.0;

    if (d < 1e-10) {
        // Objects at same location
        result.max_probability = 1.0;
        result.dilution_threshold_km = 0.0;
        result.sigma_star_km = 0.0;
        return result;
    }

    if (Rc >= d) {
        // Hard bodies overlap
        result.max_probability = 1.0;
        result.dilution_threshold_km = 0.0;
        result.sigma_star_km = 0.0;
        return result;
    }

    // Spherical case (SOCRATES default: AR = 1, fixed covariance orientation)
    // σ* = d / √2  (sigma at which probability is maximized)
    double sigma_star = d / std::sqrt(2.0);

    // P_max = (Rc² / d²) × e⁻¹
    double Pmax = (Rc * Rc) / (d * d) * std::exp(-1.0);

    // Clamp to [0, 1]
    if (Pmax > 1.0) Pmax = 1.0;

    result.max_probability = Pmax;
    result.sigma_star_km = sigma_star / 1000.0;
    result.dilution_threshold_km = sigma_star / 1000.0;

    return result;
}

// ============================================================================
// Full Collision Probability with Covariance
// ============================================================================

double collision_probability(
    const StateVector& state1, const StateVector& state2,
    const double cov1[9], const double cov2[9],
    double combined_radius_km)
{
    // Relative position and velocity
    double dr[3] = { state1.x - state2.x, state1.y - state2.y, state1.z - state2.z };
    double dv[3] = { state1.vx - state2.vx, state1.vy - state2.vy, state1.vz - state2.vz };

    double d = norm(dr);  // miss distance
    double v = norm(dv);  // relative speed

    if (v < 1e-10) return 0.0;

    // Combined covariance (uncorrelated assumption)
    double C[9];
    for (int i = 0; i < 9; i++) C[i] = cov1[i] + cov2[i];

    // Build collision plane basis (⊥ to relative velocity)
    double zhat[3] = { dv[0] / v, dv[1] / v, dv[2] / v };

    // Find a vector not parallel to zhat
    double temp[3] = { 1, 0, 0 };
    if (std::abs(zhat[0]) > 0.9) { temp[0] = 0; temp[1] = 1; }

    // xhat = normalize(temp - (temp·zhat)zhat)
    double d_tz = dot(temp, zhat);
    double xhat[3] = { temp[0] - d_tz * zhat[0], temp[1] - d_tz * zhat[1], temp[2] - d_tz * zhat[2] };
    normalize(xhat);

    // yhat = zhat × xhat
    double yhat[3];
    cross(zhat, xhat, yhat);

    // Project covariance onto collision plane (2D)
    // C_2D = P * C * P^T where P = [xhat; yhat] (2×3)
    // C_2D[0][0] = xhat^T * C * xhat
    // C_2D[0][1] = xhat^T * C * yhat
    // C_2D[1][1] = yhat^T * C * yhat

    auto matmul_vec = [](const double M[9], const double v[3], double out[3]) {
        out[0] = M[0]*v[0] + M[1]*v[1] + M[2]*v[2];
        out[1] = M[3]*v[0] + M[4]*v[1] + M[5]*v[2];
        out[2] = M[6]*v[0] + M[7]*v[1] + M[8]*v[2];
    };

    double Cx[3], Cy[3];
    matmul_vec(C, xhat, Cx);
    matmul_vec(C, yhat, Cy);

    double sxx = dot(xhat, Cx);  // σ²_xx
    double sxy = dot(xhat, Cy);  // σ²_xy
    double syy = dot(yhat, Cy);  // σ²_yy

    // Project miss vector onto collision plane
    double xm = dot(dr, xhat);
    double ym = dot(dr, yhat);
    double miss_2d = std::sqrt(xm * xm + ym * ym);

    // Eigendecompose 2×2 covariance
    double trace = sxx + syy;
    double det = sxx * syy - sxy * sxy;
    double disc = std::sqrt(std::max(0.0, trace * trace / 4.0 - det));
    double sigma1_sq = trace / 2.0 + disc;  // Larger eigenvalue
    double sigma2_sq = trace / 2.0 - disc;  // Smaller eigenvalue

    if (sigma1_sq <= 0 || sigma2_sq <= 0) {
        // Degenerate covariance — fall back to Alfano spherical
        auto result = alfano_max_probability(d, combined_radius_km);
        return result.max_probability;
    }

    double sigma1 = std::sqrt(sigma1_sq);
    double sigma2 = std::sqrt(sigma2_sq);

    // For now, use the small hard-body approximation (Foster & Estes 1992)
    // P_c ≈ (Rc² / (2σ₁σ₂)) × exp(-0.5 × (xm²/σ₁² + ym²/σ₂²))
    double Rc = combined_radius_km;
    double Pc = (Rc * Rc) / (2.0 * sigma1 * sigma2) *
                std::exp(-0.5 * (xm * xm / sigma1_sq + ym * ym / sigma2_sq));

    return std::min(Pc, 1.0);
}

// ============================================================================
// Full Conjunction Assessment
// ============================================================================

ConjunctionEvent assess_conjunction(
    const TLE& tle1, const TLE& tle2,
    double start_jd, double duration_days,
    double radius1_m, double radius2_m)
{
    ConjunctionEvent event;
    event.obj1 = tle1;
    event.obj2 = tle2;

    // Find TCA
    event.tca_jd = find_tca(tle1, tle2, start_jd, duration_days);
    event.tca_iso = jd_to_iso(event.tca_jd);

    // Propagate to TCA
    event.state1 = propagate_sgp4(tle1, event.tca_jd);
    event.state2 = propagate_sgp4(tle2, event.tca_jd);

    // Compute range and relative speed
    double dx = event.state1.x - event.state2.x;
    double dy = event.state1.y - event.state2.y;
    double dz = event.state1.z - event.state2.z;
    event.min_range_km = std::sqrt(dx * dx + dy * dy + dz * dz);

    double dvx = event.state1.vx - event.state2.vx;
    double dvy = event.state1.vy - event.state2.vy;
    double dvz = event.state1.vz - event.state2.vz;
    event.rel_speed_kms = std::sqrt(dvx * dvx + dvy * dvy + dvz * dvz);

    // RTN relative position
    inertial_to_rtn(event.state1, event.state2,
                    event.rel_pos_r, event.rel_pos_t, event.rel_pos_n,
                    event.rel_vel_r, event.rel_vel_t, event.rel_vel_n);

    // Days since epoch
    event.dse1 = event.tca_jd - tle1.epoch_jd;
    event.dse2 = event.tca_jd - tle2.epoch_jd;

    // Default covariance (SOCRATES-style, meters → km)
    event.cov_r1 = DEFAULT_COV_R_M / 1000.0;
    event.cov_t1 = DEFAULT_COV_T_M / 1000.0;
    event.cov_n1 = DEFAULT_COV_N_M / 1000.0;
    event.cov_r2 = DEFAULT_COV_R_M / 1000.0;
    event.cov_t2 = DEFAULT_COV_T_M / 1000.0;
    event.cov_n2 = DEFAULT_COV_N_M / 1000.0;

    // Alfano maximum probability
    double combined_radius_km = (radius1_m + radius2_m) / 1000.0;
    auto prob = alfano_max_probability(event.min_range_km, combined_radius_km);
    event.max_probability = prob.max_probability;
    event.dilution_threshold_km = prob.dilution_threshold_km;

    return event;
}

// ============================================================================
// Conjunction Screening
// ============================================================================

std::vector<ConjunctionEvent> screen_conjunctions(
    const std::vector<TLE>& primary_tles,
    const std::vector<TLE>& secondary_tles,
    double start_jd,
    double duration_days,
    double threshold_km)
{
    std::vector<ConjunctionEvent> events;

    for (const auto& primary : primary_tles) {
        for (const auto& secondary : secondary_tles) {
            if (primary.norad_cat_id == secondary.norad_cat_id) continue;

            try {
                auto event = assess_conjunction(primary, secondary, start_jd, duration_days);
                if (event.min_range_km <= threshold_km) {
                    events.push_back(event);
                }
            } catch (...) {
                // SGP4 error, skip this pair
            }
        }
    }

    // Sort by max probability descending
    std::sort(events.begin(), events.end(),
              [](const ConjunctionEvent& a, const ConjunctionEvent& b) {
                  return a.max_probability > b.max_probability;
              });

    return events;
}

} // namespace conjunction
