
cmake_minimum_required(VERSION 3.5...3.30)

project(CurveKID)

# --------------------------------------------------
# Geogram source directory
# --------------------------------------------------

if("$ENV{GEOGRAM_SOURCE_DIR}" STREQUAL "")
    set(
        GEOGRAM_SOURCE_DIR
        "${CMAKE_CURRENT_SOURCE_DIR}/../../../geogram/"
        CACHE PATH "Geogram source directory"
    )
else()
    set(
        GEOGRAM_SOURCE_DIR
        "$ENV{GEOGRAM_SOURCE_DIR}"
        CACHE PATH "Geogram source directory"
    )
endif()

# --------------------------------------------------
# Graphite source directory
# --------------------------------------------------

if(NOT "$ENV{GRAPHITE_SOURCE_DIR}" STREQUAL "")
    set(
        GRAPHITE_SOURCE_DIR
        "$ENV{GRAPHITE_SOURCE_DIR}"
        CACHE PATH "Graphite source directory"
    )
else()
    set(
        GRAPHITE_SOURCE_DIR
        "${CMAKE_CURRENT_SOURCE_DIR}/../../../"
        CACHE PATH "Graphite source directory"
    )
endif()

# --------------------------------------------------
# Graphite CMake configuration
# --------------------------------------------------

include("${GRAPHITE_SOURCE_DIR}/cmake/graphite.cmake")

if(EXISTS "${GRAPHITE_SOURCE_DIR}/CMakeOptions.txt")
    include("${GRAPHITE_SOURCE_DIR}/CMakeOptions.txt")
endif()

# --------------------------------------------------
# Include directories
# --------------------------------------------------

include_directories(
    "${GRAPHITE_SOURCE_DIR}/plugins"
)

# --------------------------------------------------
# Output directories
# --------------------------------------------------

set(
    CMAKE_RUNTIME_OUTPUT_DIRECTORY
    "${GRAPHITE_SOURCE_DIR}/${RELATIVE_BIN_OUTPUT_DIR}"
)

set(
    CMAKE_LIBRARY_OUTPUT_DIRECTORY
    "${GRAPHITE_SOURCE_DIR}/${RELATIVE_LIB_OUTPUT_DIR}"
)

# --------------------------------------------------
# Source files
# --------------------------------------------------

aux_source_directories(
    SOURCES
    "Source Files\\common"
    common
)

aux_source_directories(
    SOURCES
    "Source Files\\algo"
    algo
)

aux_source_directories(
    SOURCES
    "Source Files\\commands"
    commands
)

# --------------------------------------------------
# GOM code generation
# --------------------------------------------------

gomgen(CurveKID)

# --------------------------------------------------
# Build CurveKID plugin
# --------------------------------------------------

add_library(CurveKID ${SOURCES})

target_link_libraries(
    CurveKID
    mesh_gfx
    geogram
)

# --------------------------------------------------
# Visual Studio project organization
# --------------------------------------------------

set_target_properties(
    CurveKID
    PROPERTIES
    FOLDER "GRAPHITE/Plugins"
)
