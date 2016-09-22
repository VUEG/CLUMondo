cmake_minimum_required(VERSION 3.6)
project(CLUMondo)

# The version number.
set(CLUMondo_VERSION_MAJOR 1)
set(CLUMondo_VERSION_MINOR 1)
set(CLUMondo_VERSION_PATCH 0)
set(CLUMondo_VERSION "${CLUMondo_VERSION_MAJOR}.${CLUMondo_VERSION_MINOR}.${CLUMondo_VERSION_PATCH}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build/lib")

include_directories(include)
add_subdirectory(clumondo)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
        "include/CLUMondoConfig.h.in"
        "include/CLUMondoConfig.h"
)
