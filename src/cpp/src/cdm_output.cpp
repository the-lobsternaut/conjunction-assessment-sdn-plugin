/**
 * CDM FlatBuffers Output — Serialize ConjunctionEvent to CCSDS CDM
 *
 * Uses spacedatastandards.org CDM schema with $CDM file identifier.
 * Output is aligned FlatBuffers binary suitable for SDN wire format.
 */

#include "conjunction/conjunction_assessment.h"
#include "conjunction/cdm_generated.h"
#include "flatbuffers/flatbuffers.h"

#include <cstring>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <cmath>

namespace conjunction {

// jd_to_iso() is declared in sgp4_propagator.h and defined in conjunction_assessment.cpp

static std::string now_iso() {
    time_t now = time(nullptr);
    struct tm* gmt = gmtime(&now);
    char buf[64];
    snprintf(buf, sizeof(buf), "%04d-%02d-%02dT%02d:%02d:%02dZ",
             gmt->tm_year + 1900, gmt->tm_mon + 1, gmt->tm_mday,
             gmt->tm_hour, gmt->tm_min, gmt->tm_sec);
    return buf;
}

// ── Build CDMObject for one conjunction participant ──

static flatbuffers::Offset<CDMObject> build_cdm_object(
    flatbuffers::FlatBufferBuilder& builder,
    const TLE& tle,
    const StateVector& state,
    const ConjunctionEvent& event,
    int obj_num,
    double cov_r, double cov_t, double cov_n) {

    // RTN position/velocity relative to TCA
    double pos_r, pos_t, pos_n, vel_r, vel_t, vel_n;
    if (obj_num == 1) {
        pos_r = event.rel_pos_r; pos_t = event.rel_pos_t; pos_n = event.rel_pos_n;
        vel_r = event.rel_vel_r; vel_t = event.rel_vel_t; vel_n = event.rel_vel_n;
    } else {
        pos_r = 0; pos_t = 0; pos_n = 0;  // Object 2 is at origin in relative frame
        vel_r = 0; vel_t = 0; vel_n = 0;
    }

    // Build CAT (catalog) object
    auto obj_name = builder.CreateString(tle.name);
    auto obj_id = builder.CreateString(
        std::to_string(tle.norad_cat_id > 0 ? tle.norad_cat_id : 0));

    CATBuilder cat_builder(builder);
    cat_builder.add_OBJECT_NAME(obj_name);
    cat_builder.add_NORAD_CAT_ID(tle.norad_cat_id);
    auto cat = cat_builder.Finish();

    // Build covariance array (6×6 lower triangular = 21 elements for pos+vel)
    // CDM uses 9×9 (45 elements) but we only fill the 6×6 block
    std::vector<double> cov_data(45, 0.0);
    // CR_R (index 0)
    cov_data[0] = cov_r * cov_r / 1e6;  // m² → km²
    // CT_T (index 2)
    cov_data[2] = cov_t * cov_t / 1e6;
    // CN_N (index 5)
    cov_data[5] = cov_n * cov_n / 1e6;

    auto cov_vec = builder.CreateVector(cov_data);

    auto comment = builder.CreateString("SGP4 propagation");
    auto gravity = builder.CreateString("SGP4/SDP4");
    auto atm = builder.CreateString("None");

    CDMObjectBuilder obj_builder(builder);
    obj_builder.add_COMMENT(comment);
    obj_builder.add_OBJECT(cat);
    obj_builder.add_GRAVITY_MODEL(gravity);
    obj_builder.add_ATMOSPHERIC_MODEL(atm);
    obj_builder.add_X(pos_r);
    obj_builder.add_Y(pos_t);
    obj_builder.add_Z(pos_n);
    obj_builder.add_X_DOT(vel_r);
    obj_builder.add_Y_DOT(vel_t);
    obj_builder.add_Z_DOT(vel_n);
    obj_builder.add_COVARIANCE(cov_vec);

    return obj_builder.Finish();
}

// ── Public API: Serialize conjunction event to CDM FlatBuffers ──

int32_t conjunction_to_cdm(
    const ConjunctionEvent& event,
    uint8_t* output, uint32_t output_capacity) {

    flatbuffers::FlatBufferBuilder builder(4096);

    // Header strings
    auto creation_date = builder.CreateString(now_iso());
    auto originator = builder.CreateString("LOBSTERNAUT-CA");
    auto message_id = builder.CreateString(
        "CDM-" + std::to_string(event.obj1.norad_cat_id) + "-" +
        std::to_string(event.obj2.norad_cat_id));
    auto tca_str = builder.CreateString(
        event.tca_iso.empty() ? jd_to_iso(event.tca_jd) : event.tca_iso);
    auto prob_method = builder.CreateString(event.probability_method);

    // Screen period
    auto screen_start = builder.CreateString(jd_to_iso(event.tca_jd - 3.5));
    auto screen_stop = builder.CreateString(jd_to_iso(event.tca_jd + 3.5));

    // Build object entries
    auto obj1 = build_cdm_object(builder, event.obj1, event.state1, event,
                                  1, event.cov_r1, event.cov_t1, event.cov_n1);
    auto obj2 = build_cdm_object(builder, event.obj2, event.state2, event,
                                  2, event.cov_r2, event.cov_t2, event.cov_n2);

    // Build CDM
    CDMBuilder cdm_builder(builder);
    cdm_builder.add_CCSDS_CDM_VERS(1.0);
    cdm_builder.add_CREATION_DATE(creation_date);
    cdm_builder.add_ORIGINATOR(originator);
    cdm_builder.add_MESSAGE_ID(message_id);
    cdm_builder.add_TCA(tca_str);
    cdm_builder.add_MISS_DISTANCE(event.min_range_km);
    cdm_builder.add_RELATIVE_SPEED(event.rel_speed_kms);
    cdm_builder.add_RELATIVE_POSITION_R(event.rel_pos_r);
    cdm_builder.add_RELATIVE_POSITION_T(event.rel_pos_t);
    cdm_builder.add_RELATIVE_POSITION_N(event.rel_pos_n);
    cdm_builder.add_RELATIVE_VELOCITY_R(event.rel_vel_r);
    cdm_builder.add_RELATIVE_VELOCITY_T(event.rel_vel_t);
    cdm_builder.add_RELATIVE_VELOCITY_N(event.rel_vel_n);
    cdm_builder.add_START_SCREEN_PERIOD(screen_start);
    cdm_builder.add_STOP_SCREEN_PERIOD(screen_stop);
    cdm_builder.add_SCREEN_VOLUME_SHAPE(screeningVolumeShape_ELLIPSOID);
    cdm_builder.add_SCREEN_VOLUME_X(DEFAULT_THRESHOLD_KM);
    cdm_builder.add_SCREEN_VOLUME_Y(DEFAULT_THRESHOLD_KM);
    cdm_builder.add_SCREEN_VOLUME_Z(DEFAULT_THRESHOLD_KM);
    cdm_builder.add_COLLISION_PROBABILITY(event.max_probability);
    cdm_builder.add_COLLISION_PROBABILITY_METHOD(prob_method);
    cdm_builder.add_OBJECT1(obj1);
    cdm_builder.add_OBJECT2(obj2);

    auto cdm = cdm_builder.Finish();
    builder.Finish(cdm, "$CDM");

    // Copy to output
    auto buf = builder.GetBufferPointer();
    auto size = builder.GetSize();

    if (size > output_capacity) return -2;  // buffer too small
    std::memcpy(output, buf, size);
    return static_cast<int32_t>(size);
}

// ── Batch: multiple conjunctions to CDM collection ──

int32_t conjunctions_to_cdm_batch(
    const std::vector<ConjunctionEvent>& events,
    uint8_t* output, uint32_t output_capacity) {

    // For a collection, we serialize each CDM individually and pack them
    // sequentially with size prefixes (standard FlatBuffers size-prefixed format)
    uint32_t offset = 0;

    for (const auto& event : events) {
        if (offset + 4 >= output_capacity) return -2;

        // Reserve space for size prefix
        uint32_t remaining = output_capacity - offset - 4;
        int32_t written = conjunction_to_cdm(event, output + offset + 4, remaining);

        if (written < 0) return written;  // propagate error

        // Write size prefix (little-endian)
        uint32_t size = static_cast<uint32_t>(written);
        std::memcpy(output + offset, &size, 4);
        offset += 4 + size;
    }

    return static_cast<int32_t>(offset);
}

} // namespace conjunction
