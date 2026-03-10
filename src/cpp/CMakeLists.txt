cmake_minimum_required(VERSION 3.15)
project(conjunction_assessment LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# dnwrnr/sgp4 library
set(SGP4_DIR "${CMAKE_SOURCE_DIR}/deps/sgp4/libsgp4")

file(GLOB SGP4_SOURCES "${SGP4_DIR}/*.cc")
add_library(sgp4_lib STATIC ${SGP4_SOURCES})
target_include_directories(sgp4_lib PUBLIC "${SGP4_DIR}")

# SDS FlatBuffers (for WASM CDM output)
if(NOT DEFINED SDS_ROOT)
    if(DEFINED ENV{SDS_ROOT})
        set(SDS_ROOT "$ENV{SDS_ROOT}")
    else()
        foreach(_candidate
            "${CMAKE_SOURCE_DIR}/../../../../space-data-network-project"
            "$ENV{HOME}/software/space-data-network-project"
        )
            if(EXISTS "${_candidate}/packages/flatbuffers/include")
                get_filename_component(SDS_ROOT "${_candidate}" ABSOLUTE)
                break()
            endif()
        endforeach()
    endif()
endif()

if(SDS_ROOT)
    message(STATUS "SDS_ROOT: ${SDS_ROOT}")
endif()

# Conjunction assessment library
set(CA_SOURCES src/conjunction_assessment.cpp)

add_library(conjunction_lib STATIC ${CA_SOURCES})
target_include_directories(conjunction_lib PUBLIC include "${SGP4_DIR}")
target_link_libraries(conjunction_lib PUBLIC sgp4_lib)

# ============================================================================
# Native test executable
# ============================================================================

if(NOT EMSCRIPTEN)
    add_executable(test_socrates tests/test_socrates.cpp)
    target_link_libraries(test_socrates conjunction_lib)
    target_include_directories(test_socrates PRIVATE include "${SGP4_DIR}")
endif()

# ============================================================================
# WASM build (future — needs CDM FlatBuffers header)
# ============================================================================

if(EMSCRIPTEN AND SDS_ROOT)
    # TODO: Add wasm_api.cpp with CDM FlatBuffers output
    message(STATUS "WASM build: CDM output will be added later")
endif()
