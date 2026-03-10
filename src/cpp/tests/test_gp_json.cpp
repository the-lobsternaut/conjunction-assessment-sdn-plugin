/**
 * Test: GP JSON/CSV Parser + GP → TLE Conversion
 *
 * Validates:
 *   1. CelesTrak JSON GP parsing
 *   2. CelesTrak CSV GP parsing
 *   3. GP → TLE conversion produces valid TLEs for SGP4
 *   4. Propagation from GP-derived TLE matches original TLE propagation
 */

#include "conjunction/gp_json.h"
#include "conjunction/conjunction_assessment.h"
#include "conjunction/screening.h"
#include "Tle.h"
#include "SGP4.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cassert>

using namespace conjunction;

static int tests_passed = 0;
static int tests_failed = 0;

#define CHECK(cond, msg) do { \
    if (cond) { tests_passed++; std::cout << "  ✓ " << msg << std::endl; } \
    else { tests_failed++; std::cout << "  ✗ " << msg << std::endl; } \
} while(0)

void test_json_parsing() {
    std::cout << "\n=== GP JSON Parsing ===" << std::endl;

    // SOCRATES-style JSON (two objects)
    const char* json = R"([
        {"OBJECT_NAME":"STARLINK-32469","OBJECT_ID":"2024-196Q",
         "EPOCH":"2026-03-09T09:20:49.941888","MEAN_MOTION":15.30192485,
         "ECCENTRICITY":0.0001054,"INCLINATION":53.1601,
         "RA_OF_ASC_NODE":102.5175,"ARG_OF_PERICENTER":81.5233,
         "MEAN_ANOMALY":278.5887,"EPHEMERIS_TYPE":0,
         "CLASSIFICATION_TYPE":"U","NORAD_CAT_ID":61721,
         "ELEMENT_SET_NO":999,"REV_AT_EPOCH":7711,
         "BSTAR":3.1791e-5,"MEAN_MOTION_DOT":5.82e-6,"MEAN_MOTION_DDOT":0},
        {"OBJECT_NAME":"OBJECT BC","OBJECT_ID":"2025-313BC",
         "EPOCH":"2026-03-08T17:33:30.106080","MEAN_MOTION":15.29053282,
         "ECCENTRICITY":0.00085081,"INCLINATION":97.3987,
         "RA_OF_ASC_NODE":144.3551,"ARG_OF_PERICENTER":331.493,
         "MEAN_ANOMALY":28.5844,"EPHEMERIS_TYPE":0,
         "CLASSIFICATION_TYPE":"U","NORAD_CAT_ID":67298,
         "ELEMENT_SET_NO":999,"REV_AT_EPOCH":782,
         "BSTAR":0.0039787182,"MEAN_MOTION_DOT":0.00113849,"MEAN_MOTION_DDOT":0}
    ])";

    auto elements = parse_gp_json(json);
    CHECK(elements.size() == 2, "Parsed 2 GP elements from JSON");

    if (elements.size() >= 2) {
        CHECK(elements[0].object_name == "STARLINK-32469", "Object 1 name");
        CHECK(elements[0].norad_cat_id == 61721, "Object 1 NORAD ID");
        CHECK(std::abs(elements[0].mean_motion - 15.30192485) < 1e-6, "Object 1 mean motion");
        CHECK(std::abs(elements[0].inclination - 53.1601) < 1e-3, "Object 1 inclination");
        CHECK(elements[0].perigee_km > 0 && elements[0].perigee_km < 800, "Object 1 perigee altitude");
        CHECK(elements[0].apogee_km > 0 && elements[0].apogee_km < 800, "Object 1 apogee altitude");

        CHECK(elements[1].object_name == "OBJECT BC", "Object 2 name");
        CHECK(elements[1].norad_cat_id == 67298, "Object 2 NORAD ID");
        CHECK(std::abs(elements[1].eccentricity - 0.00085081) < 1e-7, "Object 2 eccentricity");
    }
}

void test_csv_parsing() {
    std::cout << "\n=== GP CSV Parsing ===" << std::endl;

    const char* csv =
        "OBJECT_NAME,OBJECT_ID,EPOCH,MEAN_MOTION,ECCENTRICITY,INCLINATION,"
        "RA_OF_ASC_NODE,ARG_OF_PERICENTER,MEAN_ANOMALY,EPHEMERIS_TYPE,"
        "CLASSIFICATION_TYPE,NORAD_CAT_ID,ELEMENT_SET_NO,REV_AT_EPOCH,"
        "BSTAR,MEAN_MOTION_DOT,MEAN_MOTION_DDOT\n"
        "ISS (ZARYA),1998-067A,2026-03-10T03:48:27.801504,15.48557200,"
        ".00079716,51.6324,70.0375,176.2903,183.8146,0,U,25544,999,55640,"
        ".17625562E-3,.913E-4,0\n";

    auto elements = parse_gp_csv(csv);
    CHECK(elements.size() == 1, "Parsed 1 GP element from CSV");

    if (!elements.empty()) {
        CHECK(elements[0].object_name == "ISS (ZARYA)", "ISS name");
        CHECK(elements[0].norad_cat_id == 25544, "ISS NORAD ID");
        CHECK(std::abs(elements[0].mean_motion - 15.485572) < 1e-4, "ISS mean motion");
        CHECK(elements[0].perigee_km > 350 && elements[0].perigee_km < 450, "ISS perigee ~410km");
        CHECK(elements[0].apogee_km > 350 && elements[0].apogee_km < 450, "ISS apogee ~420km");
    }
}

void test_gp_to_tle_conversion() {
    std::cout << "\n=== GP → TLE Conversion ===" << std::endl;

    // Parse GP JSON
    const char* json = R"([
        {"OBJECT_NAME":"STARLINK-32469","OBJECT_ID":"2024-196Q",
         "EPOCH":"2026-03-09T09:20:49.941888","MEAN_MOTION":15.30192485,
         "ECCENTRICITY":0.0001054,"INCLINATION":53.1601,
         "RA_OF_ASC_NODE":102.5175,"ARG_OF_PERICENTER":81.5233,
         "MEAN_ANOMALY":278.5887,"EPHEMERIS_TYPE":0,
         "CLASSIFICATION_TYPE":"U","NORAD_CAT_ID":61721,
         "ELEMENT_SET_NO":999,"REV_AT_EPOCH":7711,
         "BSTAR":3.1791e-5,"MEAN_MOTION_DOT":5.82e-6,"MEAN_MOTION_DDOT":0}
    ])";

    auto elements = parse_gp_json(json);
    CHECK(elements.size() == 1, "Parsed GP element");

    if (!elements.empty()) {
        auto tle = gp_to_tle(elements[0]);
        CHECK(tle.norad_cat_id == 61721, "TLE NORAD ID matches");
        CHECK(tle.line1.length() == 69, "Line 1 is 69 chars (got " + std::to_string(tle.line1.length()) + ")");
        CHECK(tle.line2.length() == 69, "Line 2 is 69 chars (got " + std::to_string(tle.line2.length()) + ")");
        CHECK(tle.line1[0] == '1', "Line 1 starts with '1'");
        CHECK(tle.line2[0] == '2', "Line 2 starts with '2'");

        std::cout << "  Line1: [" << tle.line1 << "]" << std::endl;
        std::cout << "  Line2: [" << tle.line2 << "]" << std::endl;

        // Try to parse with dnwrnr/sgp4
        try {
            libsgp4::Tle sgp4_tle(tle.name, tle.line1, tle.line2);
            CHECK(true, "dnwrnr/sgp4 accepts the generated TLE");

            // Propagate forward 1 day
            libsgp4::SGP4 sgp4(sgp4_tle);
            auto epoch = sgp4_tle.Epoch();
            libsgp4::TimeSpan one_day(1, 0, 0, 0);
            auto target = epoch + one_day;
            auto eci = sgp4.FindPosition(target);

            double r = std::sqrt(eci.Position().x * eci.Position().x +
                                eci.Position().y * eci.Position().y +
                                eci.Position().z * eci.Position().z);
            CHECK(r > 6400 && r < 7200, "Position radius in LEO range (r=" +
                  std::to_string(r) + " km)");

        } catch (const std::exception& e) {
            CHECK(false, "dnwrnr/sgp4 parse: " + std::string(e.what()));
        }
    }
}

void test_gp_json_propagation_accuracy() {
    std::cout << "\n=== GP JSON Propagation Accuracy ===" << std::endl;

    // Parse the SOCRATES pair from JSON and compare with TLE-based propagation
    const char* json = R"([
        {"OBJECT_NAME":"STARLINK-32469","OBJECT_ID":"2024-196Q",
         "EPOCH":"2026-03-09T09:20:49.941888","MEAN_MOTION":15.30192485,
         "ECCENTRICITY":0.0001054,"INCLINATION":53.1601,
         "RA_OF_ASC_NODE":102.5175,"ARG_OF_PERICENTER":81.5233,
         "MEAN_ANOMALY":278.5887,"EPHEMERIS_TYPE":0,
         "CLASSIFICATION_TYPE":"U","NORAD_CAT_ID":61721,
         "ELEMENT_SET_NO":999,"REV_AT_EPOCH":7711,
         "BSTAR":3.1791e-5,"MEAN_MOTION_DOT":5.82e-6,"MEAN_MOTION_DDOT":0},
        {"OBJECT_NAME":"OBJECT BC","OBJECT_ID":"2025-313BC",
         "EPOCH":"2026-03-08T17:33:30.106080","MEAN_MOTION":15.29053282,
         "ECCENTRICITY":0.00085081,"INCLINATION":97.3987,
         "RA_OF_ASC_NODE":144.3551,"ARG_OF_PERICENTER":331.493,
         "MEAN_ANOMALY":28.5844,"EPHEMERIS_TYPE":0,
         "CLASSIFICATION_TYPE":"U","NORAD_CAT_ID":67298,
         "ELEMENT_SET_NO":999,"REV_AT_EPOCH":782,
         "BSTAR":0.0039787182,"MEAN_MOTION_DOT":0.00113849,"MEAN_MOTION_DDOT":0}
    ])";

    auto elements = parse_gp_json(json);
    CHECK(elements.size() == 2, "Parsed 2 GP elements");

    if (elements.size() >= 2) {
        auto tle1 = gp_to_tle(elements[0]);
        auto tle2 = gp_to_tle(elements[1]);

        // Run conjunction assessment
        double start_jd = iso_to_jd("2026-03-09T18:00:00Z");
        try {
            auto event = assess_conjunction(tle1, tle2, start_jd, 7.0);
            std::cout << "  TCA: " << event.tca_iso << std::endl;
            std::cout << "  Range: " << event.min_range_km << " km" << std::endl;
            std::cout << "  Speed: " << event.rel_speed_kms << " km/s" << std::endl;
            std::cout << "  MaxProb: " << event.max_probability << std::endl;

            // SOCRATES reference: TCA=2026-03-14T11:08:50.281Z, range=0.010, speed=9.757
            double ref_tca_jd = iso_to_jd("2026-03-14T11:08:50.281Z");
            double tca_diff = std::abs(event.tca_jd - ref_tca_jd) * 86400.0;

            CHECK(tca_diff < 120.0, "TCA within 2 minutes of SOCRATES (Δ=" +
                  std::to_string(tca_diff) + "s)");
            CHECK(event.min_range_km < 0.5, "Range within 500m of SOCRATES");
            CHECK(std::abs(event.rel_speed_kms - 9.757) < 1.0, "Speed within 1 km/s of SOCRATES");
        } catch (const std::exception& e) {
            CHECK(false, "Conjunction assessment: " + std::string(e.what()));
        }
    }
}

void test_altitude_prefilter() {
    std::cout << "\n=== Altitude Prefilter ===" << std::endl;

    GPElement leo1, leo2, geo;

    // Two LEO objects at ~550 km
    leo1.perigee_km = 540; leo1.apogee_km = 560;
    leo2.perigee_km = 545; leo2.apogee_km = 555;

    // GEO object at ~35,786 km
    geo.perigee_km = 35780; geo.apogee_km = 35790;

    CHECK(altitude_overlap(leo1, leo2), "LEO1 and LEO2 overlap");
    CHECK(!altitude_overlap(leo1, geo), "LEO1 and GEO do NOT overlap");
    CHECK(!altitude_overlap(leo2, geo), "LEO2 and GEO do NOT overlap");

    // HEO (Molniya-like): 500 km perigee, 40,000 km apogee
    GPElement heo;
    heo.perigee_km = 500; heo.apogee_km = 40000;
    CHECK(altitude_overlap(leo1, heo), "LEO and HEO overlap (HEO passes through LEO)");
    CHECK(altitude_overlap(geo, heo), "GEO and HEO overlap (HEO passes through GEO)");
}

void test_alpha5_encoding() {
    std::cout << "\n--- Test: Alpha-5 NORAD ID Encoding ---" << std::endl;

    // Standard 5-digit IDs
    {
        GPElement gp{};
        gp.norad_cat_id = 25544;
        gp.classification_type = 'U';
        gp.object_id = "1998-067A";
        gp.epoch_iso = "2026-03-09T12:00:00.000000";
        gp.mean_motion = 15.5;
        gp.eccentricity = 0.0001;
        gp.inclination = 51.6;
        gp.ra_of_asc_node = 100.0;
        gp.arg_of_pericenter = 50.0;
        gp.mean_anomaly = 200.0;
        gp.bstar = 0.00005;
        gp.mean_motion_dot = 0.00001;
        gp.mean_motion_ddot = 0.0;
        gp.ephemeris_type = 0;
        gp.element_set_no = 999;
        gp.rev_at_epoch = 1000;
        gp.epoch_jd = iso_to_jd(gp.epoch_iso);

        TLE tle = gp_to_tle(gp);
        CHECK(tle.line1.substr(2, 5) == "25544", "Standard ID 25544 in line1");
        CHECK(tle.line2.substr(2, 5) == "25544", "Standard ID 25544 in line2");

        // Verify SGP4 can parse it
        try {
            libsgp4::Tle sgp4_tle("ISS", tle.line1, tle.line2);
            CHECK(sgp4_tle.NoradNumber() == 25544, "SGP4 parsed NORAD 25544");
        } catch (const std::exception& e) {
            CHECK(false, std::string("SGP4 parse standard ID: ") + e.what());
        }
    }

    // Alpha-5: 100000 → A0000
    {
        GPElement gp{};
        gp.norad_cat_id = 100000;
        gp.classification_type = 'U';
        gp.object_id = "2025-001A";
        gp.epoch_iso = "2026-03-09T12:00:00.000000";
        gp.mean_motion = 15.0;
        gp.eccentricity = 0.001;
        gp.inclination = 97.0;
        gp.ra_of_asc_node = 150.0;
        gp.arg_of_pericenter = 60.0;
        gp.mean_anomaly = 300.0;
        gp.bstar = 0.0001;
        gp.mean_motion_dot = 0.0;
        gp.mean_motion_ddot = 0.0;
        gp.ephemeris_type = 0;
        gp.element_set_no = 999;
        gp.rev_at_epoch = 100;
        gp.epoch_jd = iso_to_jd(gp.epoch_iso);

        TLE tle = gp_to_tle(gp);
        CHECK(tle.line1.substr(2, 5) == "A0000", "Alpha-5 ID 100000 → A0000 in line1");
        CHECK(tle.line2.substr(2, 5) == "A0000", "Alpha-5 ID 100000 → A0000 in line2");

        try {
            libsgp4::Tle sgp4_tle("SAT", tle.line1, tle.line2);
            CHECK(sgp4_tle.NoradNumber() == 100000, "SGP4 parsed Alpha-5 NORAD 100000");
        } catch (const std::exception& e) {
            CHECK(false, std::string("SGP4 parse Alpha-5 100000: ") + e.what());
        }
    }

    // Alpha-5: 270086 → R0086 (27 → R: P=23, Q=24, R=25... wait 27-23=4, P+4=T? No...)
    // 270086: first = 27, rest = 0086
    // 27 >= 23 → 'P' + (27-23) = 'T'
    // So 270086 → T0086
    {
        GPElement gp{};
        gp.norad_cat_id = 270086;
        gp.classification_type = 'U';
        gp.object_id = "2025-999A";
        gp.epoch_iso = "2026-03-09T12:00:00.000000";
        gp.mean_motion = 14.5;
        gp.eccentricity = 0.002;
        gp.inclination = 45.0;
        gp.ra_of_asc_node = 200.0;
        gp.arg_of_pericenter = 90.0;
        gp.mean_anomaly = 180.0;
        gp.bstar = 0.0002;
        gp.mean_motion_dot = 0.0;
        gp.mean_motion_ddot = 0.0;
        gp.ephemeris_type = 0;
        gp.element_set_no = 999;
        gp.rev_at_epoch = 50;
        gp.epoch_jd = iso_to_jd(gp.epoch_iso);

        TLE tle = gp_to_tle(gp);
        CHECK(tle.line1.substr(2, 5) == "T0086", "Alpha-5 ID 270086 → T0086 in line1");
        CHECK(tle.line2.substr(2, 5) == "T0086", "Alpha-5 ID 270086 → T0086 in line2");

        try {
            libsgp4::Tle sgp4_tle("SAT", tle.line1, tle.line2);
            CHECK(sgp4_tle.NoradNumber() == 270086, "SGP4 parsed Alpha-5 NORAD 270086");
        } catch (const std::exception& e) {
            CHECK(false, std::string("SGP4 parse Alpha-5 270086: ") + e.what());
        }
    }

    // Edge cases
    {
        GPElement gp{};
        gp.classification_type = 'U';
        gp.object_id = "2025-001A";
        gp.epoch_iso = "2026-03-09T12:00:00.000000";
        gp.mean_motion = 15.0; gp.eccentricity = 0.001;
        gp.inclination = 97.0; gp.ra_of_asc_node = 150.0;
        gp.arg_of_pericenter = 60.0; gp.mean_anomaly = 300.0;
        gp.bstar = 0.0001; gp.mean_motion_dot = 0.0;
        gp.mean_motion_ddot = 0.0; gp.ephemeris_type = 0;
        gp.element_set_no = 999; gp.rev_at_epoch = 100;
        gp.epoch_jd = iso_to_jd(gp.epoch_iso);

        // 339999 → Z9999 (max Alpha-5)
        gp.norad_cat_id = 339999;
        TLE tle = gp_to_tle(gp);
        CHECK(tle.line1.substr(2, 5) == "Z9999", "Alpha-5 max: 339999 → Z9999");

        try {
            libsgp4::Tle sgp4_tle("SAT", tle.line1, tle.line2);
            CHECK(sgp4_tle.NoradNumber() == 339999, "SGP4 parsed Alpha-5 max 339999");
        } catch (const std::exception& e) {
            CHECK(false, std::string("SGP4 parse Alpha-5 339999: ") + e.what());
        }

        // 180000 → J0000 (J=18, 180000 % 10000 = 0)
        gp.norad_cat_id = 180000;
        tle = gp_to_tle(gp);
        CHECK(tle.line1.substr(2, 5) == "J0000", "Alpha-5: 180000 → J0000");

        try {
            libsgp4::Tle sgp4_tle("SAT", tle.line1, tle.line2);
            CHECK(sgp4_tle.NoradNumber() == 180000, "SGP4 parsed Alpha-5 180000");
        } catch (const std::exception& e) {
            CHECK(false, std::string("SGP4 parse Alpha-5 180000: ") + e.what());
        }
    }
}

int main() {
    test_json_parsing();
    test_csv_parsing();
    test_gp_to_tle_conversion();
    test_gp_json_propagation_accuracy();
    test_altitude_prefilter();
    test_alpha5_encoding();

    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Passed: " << tests_passed << std::endl;
    std::cout << "Failed: " << tests_failed << std::endl;

    return tests_failed > 0 ? 1 : 0;
}
