cmake_minimum_required(VERSION 3.24)
include(CheckCXXCompilerFlag)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_program(CCACHE_PROGRAM ccache)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

if(CCACHE_PROGRAM)
  message(STATUS "Using compiler cache")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
endif()
project(KaHIP C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED OFF)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) 

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# if no build mode is specified build in release mode
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# tweak compiler flags
CHECK_CXX_COMPILER_FLAG(-funroll-loops COMPILER_SUPPORTS_FUNROLL_LOOPS)
if(COMPILER_SUPPORTS_FUNROLL_LOOPS)
  add_definitions(-funroll-loops)
endif()
CHECK_CXX_COMPILER_FLAG(-fno-stack-limit COMPILER_SUPPORTS_FNOSTACKLIMITS)
if(COMPILER_SUPPORTS_FNOSTACKLIMITS)
  add_definitions(-fno-stack-limit)
endif()
CHECK_CXX_COMPILER_FLAG(-Wall COMPILER_SUPPORTS_WALL)
if(COMPILER_SUPPORTS_WALL)
  add_definitions(-Wall)
endif()
CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
  add_definitions(-march=native)
endif()
CHECK_CXX_COMPILER_FLAG(-fpermissive COMPILER_SUPPORTS_FPERMISSIVE)
if(COMPILER_SUPPORTS_FPERMISSIVE)
  add_definitions(-fpermissive)
endif()

# check dependencies
find_package(MPI REQUIRED)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP support detected")
  add_definitions(${OpenMP_CXX_FLAGS})
else()
  message(WARNING "OpenMP not available, activating workaround")
  add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
  set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS "")
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/misc)
endif()

# 64 Bit option
#option(64BITMODE "64 bit mode" OFF)
option(64BITMODE "64 bit mode" ON)
if(64BITMODE)
  add_definitions("-DMODE64BITEDGES")
  add_definitions("-DPOINTER64=1")
endif()

option(64BITVERTEXMODE "64 bit mode" ON)
if(64BITVERTEXMODE)
  add_definitions("-DMODE64BITVERTICES")
endif()

# optimized output
option(OPTIMIZED_OUTPUT "optimized output" OFF)
if(OPTIMIZED_OUTPUT)
  add_definitions("-DKAFFPAOUTPUT")
endif()

include(FetchContent)

FetchContent_Declare(
        cpi
        GIT_REPOSITORY https://github.com/kurpicz/cpi.git
        GIT_TAG "origin/main"
)

FetchContent_MakeAvailable(cpi)

#include_directories(${cpi_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/app)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/argtable3-3.0.3)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/io)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/partition)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/tools)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/data_structure/compression_vectors)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extern/flatbuffers
        ${CMAKE_CURRENT_BINARY_DIR}/flatbuffers-build
        EXCLUDE_FROM_ALL)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extern/stxxl)

set(KAGEN_BUILD_APPS OFF)
set(KAGEN_BUILD_EXAMPLES OFF)
set(KAGEN_BUILD_TESTS OFF)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extern/KaGen)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS} -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc")

# add STXXL includes path
include_directories(${STXXL_INCLUDE_DIRS})

set(LIBKAFFPA_SOURCE_FILES
  lib/tools/quality_metrics.cpp
  lib/tools/random_functions.cpp)
add_library(libkaffpa OBJECT ${LIBKAFFPA_SOURCE_FILES})
TARGET_LINK_LIBRARIES(libkaffpa PRIVATE cpi)

set(ARGTABLE_SOURCE_FILES
extern/argtable3-3.0.3/argtable3.c)
add_library(argtab OBJECT ${ARGTABLE_SOURCE_FILES})

set(GRAPH_IO_SOURCE_FILES
        lib/io/graph_io_stream.h
        lib/io/graph_io_stream.cpp)
add_library(graphio OBJECT ${GRAPH_IO_SOURCE_FILES})
target_link_libraries(graphio PRIVATE cpi)

set(LIBSTREAMPARTITIONING_SOURCE_FILES
  lib/partition/onepass_partitioning/vertex_partitioning.cpp
  lib/partition/onepass_partitioning/floating_block.cpp
  lib/partition/onepass_partitioning/ldg.cpp
  lib/partition/onepass_partitioning/fennel_approx_sqrt.cpp
  lib/partition/onepass_partitioning/fennel.cpp)
add_library(libstreampartitioning OBJECT ${LIBSTREAMPARTITIONING_SOURCE_FILES})
target_link_libraries(libstreampartitioning PRIVATE cpi)

# generate targets for each binary
#add_executable(freight_graphs app/streammultisection.cpp $<TARGET_OBJECTS:libpart> $<TARGET_OBJECTS:argtab> $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping> $<TARGET_OBJECTS:libstreampartitioning> )
add_executable(stream_cpi app/streamcpi.cpp)
target_compile_definitions(stream_cpi PRIVATE "-DMODE_FREIGHT_GRAPHS")
target_link_libraries(stream_cpi
        PRIVATE cpi
        graphio
        libkaffpa
        argtab
        libstreampartitioning
        flatbuffers
        ${STXXL_LIBRARIES}
        ${OpenMP_CXX_LIBRARIES})
install(TARGETS stream_cpi DESTINATION bin)

add_executable(stream_cpi_generated app/generatedgraphpartition.cpp)
target_compile_definitions(stream_cpi_generated PRIVATE "-DMODE_FREIGHT_GRAPHS")
target_link_libraries(stream_cpi_generated
        PRIVATE cpi
        graphio
        libkaffpa
        argtab
        libstreampartitioning
        flatbuffers
        # Linking KaGen increases memory consumption by 3k bytes
        PUBLIC KaGen::KaGen
        ${STXXL_LIBRARIES}
        ${OpenMP_CXX_LIBRARIES})
install(TARGETS stream_cpi_generated DESTINATION bin)
