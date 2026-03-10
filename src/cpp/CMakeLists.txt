cmake_minimum_required(VERSION 3.15)
project(conjunction_assessment LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ============================================================================
# Threading support
# ============================================================================
find_package(Threads REQUIRED)

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

# ============================================================================
# Conjunction assessment library
# ============================================================================
set(CA_SOURCES
    src/conjunction_assessment.cpp
    src/gp_json.cpp
    src/kdtree.cpp
    src/screening.cpp
)

add_library(conjunction_lib STATIC ${CA_SOURCES})
target_include_directories(conjunction_lib PUBLIC include "${SGP4_DIR}")
target_link_libraries(conjunction_lib PUBLIC sgp4_lib Threads::Threads)

# ============================================================================
# Native test executables
# ============================================================================

if(NOT EMSCRIPTEN)
    # SOCRATES validation test
    add_executable(test_socrates tests/test_socrates.cpp)
    target_link_libraries(test_socrates conjunction_lib)
    target_include_directories(test_socrates PRIVATE include "${SGP4_DIR}")

    # GP JSON parser test
    add_executable(test_gp_json tests/test_gp_json.cpp)
    target_link_libraries(test_gp_json conjunction_lib)
    target_include_directories(test_gp_json PRIVATE include "${SGP4_DIR}")

    # KD-tree screening benchmark
    add_executable(test_screening tests/test_screening.cpp)
    target_link_libraries(test_screening conjunction_lib)
    target_include_directories(test_screening PRIVATE include "${SGP4_DIR}")
endif()

# ============================================================================
# WASM build (Emscripten + pthreads + SharedArrayBuffer)
# ============================================================================

if(EMSCRIPTEN)
    set(WASM_FLAGS
        -s EXPORTED_FUNCTIONS=['_parse','_convert','_screen','_malloc','_free']
        -s EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString']
        -s ALLOW_MEMORY_GROWTH=1
        -s INITIAL_MEMORY=67108864      # 64MB
        -s MAXIMUM_MEMORY=536870912     # 512MB
        -s USE_PTHREADS=1
        -s PTHREAD_POOL_SIZE=navigator.hardwareConcurrency
        -s SHARED_MEMORY=1
        -s MODULARIZE=1
        -s EXPORT_ES6=1
        -s EXPORT_NAME=ConjunctionWASM
        -s NO_EXIT_RUNTIME=1
        --no-entry
    )

    if(SDS_ROOT)
        # WASM API with CDM FlatBuffers output
        add_executable(conjunction_assessment_wasm
            src/wasm_api.cpp
            ${CA_SOURCES}
            ${SGP4_SOURCES}
        )
        target_include_directories(conjunction_assessment_wasm PRIVATE
            include "${SGP4_DIR}"
            "${SDS_ROOT}/packages/flatbuffers/include"
            "${SDS_ROOT}/packages/spacedatastandards.org/lib/cpp"
        )
        target_compile_options(conjunction_assessment_wasm PRIVATE -O2 -pthread)
        target_link_options(conjunction_assessment_wasm PRIVATE ${WASM_FLAGS})
        set_target_properties(conjunction_assessment_wasm PROPERTIES
            SUFFIX ".mjs"
            OUTPUT_NAME "conjunction_assessment"
        )
    endif()
endif()
