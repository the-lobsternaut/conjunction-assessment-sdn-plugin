/**
 * CDM FlatBuffers Output Test
 *
 * Verifies conjunction events serialize to valid CDM FlatBuffers
 * with $CDM file identifier.
 */

#include "conjunction/conjunction_assessment.h"
#include "conjunction/cdm_generated.h"
#include "flatbuffers/flatbuffers.h"

#include <iostream>
#include <cmath>
#include <cstring>

using namespace conjunction;

static int tests_passed = 0;
static int tests_failed = 0;

#define CHECK(cond, msg) do { \
    if (cond) { tests_passed++; std::cout << "  ✓ " << msg << std::endl; } \
    else { tests_failed++; std::cout << "  ✗ " << msg << std::endl; } \
} while(0)

void test_single_cdm_output() {
    std::cout << "\n--- Test: Single CDM Output ---" << std::endl;

    // Create a synthetic conjunction event
    ConjunctionEvent event;
    event.obj1.name = "STARLINK-1234";
    event.obj1.norad_cat_id = 45678;
    event.obj2.name = "COSMOS 2251 DEB";
    event.obj2.norad_cat_id = 34567;

    event.tca_jd = 2460784.5;
    event.tca_iso = "2025-03-15T12:00:00.000Z";
    event.min_range_km = 0.234;
    event.rel_speed_kms = 14.7;

    event.rel_pos_r = 0.1;
    event.rel_pos_t = 0.15;
    event.rel_pos_n = -0.08;
    event.rel_vel_r = 2.3;
    event.rel_vel_t = 12.1;
    event.rel_vel_n = -5.4;

    event.max_probability = 1.23e-4;
    event.probability_method = "ALFANO-MAXPROB";

    event.cov_r1 = DEFAULT_COV_R_M;
    event.cov_t1 = DEFAULT_COV_T_M;
    event.cov_n1 = DEFAULT_COV_N_M;
    event.cov_r2 = DEFAULT_COV_R_M;
    event.cov_t2 = DEFAULT_COV_T_M;
    event.cov_n2 = DEFAULT_COV_N_M;

    // Serialize
    uint8_t buffer[16384];
    int32_t written = conjunction_to_cdm(event, buffer, sizeof(buffer));

    CHECK(written > 0, "CDM serialization succeeded (" + std::to_string(written) + " bytes)");

    // Verify file identifier
    bool has_id = flatbuffers::BufferHasIdentifier(buffer, "$CDM");
    CHECK(has_id, "Buffer has $CDM file identifier");

    // Deserialize and verify fields
    auto cdm = GetCDM(buffer);
    CHECK(cdm != nullptr, "CDM deserialized successfully");

    if (cdm) {
        CHECK(cdm->CCSDS_CDM_VERS() == 1.0, "CCSDS version = 1.0");
        CHECK(std::string(cdm->ORIGINATOR()->c_str()) == "LOBSTERNAUT-CA", "Originator correct");

        CHECK(cdm->TCA() != nullptr, "TCA string present");
        if (cdm->TCA()) {
            CHECK(std::string(cdm->TCA()->c_str()) == "2025-03-15T12:00:00.000Z", "TCA matches");
        }

        CHECK(std::abs(cdm->MISS_DISTANCE() - 0.234) < 0.001, "Miss distance = 0.234 km");
        CHECK(std::abs(cdm->RELATIVE_SPEED() - 14.7) < 0.01, "Relative speed = 14.7 km/s");
        CHECK(std::abs(cdm->COLLISION_PROBABILITY() - 1.23e-4) < 1e-8, "Probability = 1.23e-4");

        // Check RTN components
        CHECK(std::abs(cdm->RELATIVE_POSITION_R() - 0.1) < 0.001, "Rel pos R correct");
        CHECK(std::abs(cdm->RELATIVE_POSITION_T() - 0.15) < 0.001, "Rel pos T correct");
        CHECK(std::abs(cdm->RELATIVE_POSITION_N() - (-0.08)) < 0.001, "Rel pos N correct");

        // Check objects
        CHECK(cdm->OBJECT1() != nullptr, "Object 1 present");
        CHECK(cdm->OBJECT2() != nullptr, "Object 2 present");

        if (cdm->OBJECT1() && cdm->OBJECT1()->OBJECT()) {
            CHECK(cdm->OBJECT1()->OBJECT()->NORAD_CAT_ID() == 45678, "Object 1 NORAD ID");
        }

        if (cdm->OBJECT2() && cdm->OBJECT2()->OBJECT()) {
            CHECK(cdm->OBJECT2()->OBJECT()->NORAD_CAT_ID() == 34567, "Object 2 NORAD ID");
        }

        // Check screening volume
        CHECK(cdm->SCREEN_VOLUME_SHAPE() == screeningVolumeShape_ELLIPSOID, "Screening volume shape = ELLIPSOID");
        CHECK(std::abs(cdm->SCREEN_VOLUME_X() - 5.0) < 0.01, "Screening volume X = 5 km");
    }

    std::cout << "    CDM size: " << written << " bytes" << std::endl;
}

void test_batch_cdm_output() {
    std::cout << "\n--- Test: Batch CDM Output ---" << std::endl;

    std::vector<ConjunctionEvent> events;
    for (int i = 0; i < 5; i++) {
        ConjunctionEvent event;
        event.obj1.name = "SAT-" + std::to_string(i);
        event.obj1.norad_cat_id = 10000 + i;
        event.obj2.name = "DEB-" + std::to_string(i);
        event.obj2.norad_cat_id = 20000 + i;
        event.tca_jd = 2460784.5 + i * 0.1;
        event.min_range_km = 1.0 + i * 0.5;
        event.rel_speed_kms = 10.0 + i;
        event.max_probability = 1e-5 * (i + 1);
        event.probability_method = "ALFANO-MAXPROB";
        event.cov_r1 = event.cov_r2 = DEFAULT_COV_R_M;
        event.cov_t1 = event.cov_t2 = DEFAULT_COV_T_M;
        event.cov_n1 = event.cov_n2 = DEFAULT_COV_N_M;
        events.push_back(event);
    }

    uint8_t buffer[65536];
    int32_t written = conjunctions_to_cdm_batch(events, buffer, sizeof(buffer));

    CHECK(written > 0, "Batch CDM serialization succeeded (" + std::to_string(written) + " bytes)");

    // Parse the size-prefixed buffer
    int count = 0;
    uint32_t offset = 0;
    while (offset < static_cast<uint32_t>(written)) {
        uint32_t size;
        std::memcpy(&size, buffer + offset, 4);
        offset += 4;

        bool valid = flatbuffers::BufferHasIdentifier(buffer + offset, "$CDM");
        CHECK(valid, "CDM #" + std::to_string(count) + " has valid $CDM identifier");

        auto cdm = GetCDM(buffer + offset);
        if (cdm) {
            CHECK(cdm->OBJECT1() != nullptr, "CDM #" + std::to_string(count) + " has Object 1");
        }

        offset += size;
        count++;
    }

    CHECK(count == 5, "Parsed 5 CDMs from batch (" + std::to_string(count) + " found)");
    std::cout << "    Batch size: " << written << " bytes for " << count << " CDMs" << std::endl;
}

void test_buffer_too_small() {
    std::cout << "\n--- Test: Buffer Too Small ---" << std::endl;

    ConjunctionEvent event;
    event.obj1.name = "SAT-A";
    event.obj1.norad_cat_id = 11111;
    event.obj2.name = "SAT-B";
    event.obj2.norad_cat_id = 22222;
    event.tca_jd = 2460784.5;
    event.min_range_km = 1.0;
    event.rel_speed_kms = 10.0;
    event.max_probability = 1e-5;
    event.probability_method = "ALFANO-MAXPROB";

    uint8_t tiny_buffer[10];
    int32_t written = conjunction_to_cdm(event, tiny_buffer, sizeof(tiny_buffer));
    CHECK(written == -2, "Returns -2 for buffer too small");
}

int main() {
    std::cout << "=== CDM FlatBuffers Output Tests ===" << std::endl;

    test_single_cdm_output();
    test_batch_cdm_output();
    test_buffer_too_small();

    std::cout << "\n=== Summary: " << tests_passed << " passed, "
              << tests_failed << " failed ===" << std::endl;

    return tests_failed > 0 ? 1 : 0;
}
