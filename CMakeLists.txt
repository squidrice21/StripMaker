cmake_minimum_required(VERSION 3.4...3.18)

project(stripmaker)

## OPTIONS -- use tbb static or dynamic
option(USE_TBB_STATIC    "Build with static TBB"   ON)

# Set c++11
# https://stackoverflow.com/questions/10851247/how-to-activate-c-11-in-cmake
if (POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW)
endif ()
set (CMAKE_CXX_STANDARD 20)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

## Add the cmake/ folder to the module path so we can use the files in there.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
# For FETCH CONTEXTS, put the dependencies in root folder to allow sharing them
# among debug/release builds
# If a value is specified when running cmake, it will override this.
set(FETCHCONTENT_BASE_DIR "${CMAKE_CURRENT_LIST_DIR}/_external" CACHE STRING "")

# Visual Studio compiler with static runtime libraries
if(MSVC AND MT)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")
endif()
if(MSVC)
    # cmake_policy(SET CMP0091 OLD)
    # set(CMAKE_MSVC_RUNTIME_LIBRARY "")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /bigobj")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /bigobj")
    set(CMAKE_ASM_MASM_COMPILE_OPTIONS_MSVC_RUNTIME_LIBRARY_MultiThreaded         "")
    set(CMAKE_ASM_MASM_COMPILE_OPTIONS_MSVC_RUNTIME_LIBRARY_MultiThreadedDLL      "")
    set(CMAKE_ASM_MASM_COMPILE_OPTIONS_MSVC_RUNTIME_LIBRARY_MultiThreadedDebug    "")
    set(CMAKE_ASM_MASM_COMPILE_OPTIONS_MSVC_RUNTIME_LIBRARY_MultiThreadedDebugDLL "")
endif()

# Source files
file(GLOB STROKE_STRIP_SOURCE_FILES stroke_strip_src/*.cpp stroke_strip_src/*.h)
list(FILTER STROKE_STRIP_SOURCE_FILES EXCLUDE REGEX ".*Main.cpp$")
# add_library(lib${CMAKE_PROJECT_NAME} ${STROKE_STRIP_SOURCE_FILES})

file(GLOB SOURCE_FILES
    src/*.cpp src/*.h
    src/feature/*.cpp src/feature/*.h
    src/sample/*.cpp src/sample/*.h
    src/measure/*.cpp src/measure/*.h
    src/classifier/*.cpp src/classifier/*.h
    src/uncutting/*.cpp src/uncutting/*.h
    src/solve/*.cpp src/solve/*.h
    src/endpoint/*.cpp src/endpoint/*.h)
list(FILTER SOURCE_FILES EXCLUDE REGEX ".*Main.*.cpp$")
list(FILTER SOURCE_FILES EXCLUDE REGEX ".*Main.*.h$")
add_library(lib${CMAKE_PROJECT_NAME} ${SOURCE_FILES} ${STROKE_STRIP_SOURCE_FILES})

target_compile_definitions(lib${CMAKE_PROJECT_NAME} PUBLIC _USE_MATH_DEFINES)

# Entrances
add_executable(clustering-preprocess "src/MainPreprocess.cpp")
add_executable(clustering-feature "src/MainFeature.cpp")
add_executable(clustering-secondary-feature "src/MainSecondaryFeature.cpp")
add_executable(clustering-sample "src/MainSample.cpp")
add_executable(clustering-secondary-sample "src/MainSecondarySample.cpp")
add_executable(clustering-solve "src/MainSolve.cpp")
add_executable(clustering-uncut "src/MainUncutting.cpp")
add_executable(clustering-endend "src/MainEndend.cpp")
add_executable(strokestrip "src/MainStrokeStrip.cpp")

# Data generation entrance
include(cli11)
target_link_libraries(clustering-preprocess PUBLIC lib${CMAKE_PROJECT_NAME} CLI11::CLI11)
target_link_libraries(clustering-feature PUBLIC lib${CMAKE_PROJECT_NAME} CLI11::CLI11)
target_link_libraries(clustering-secondary-feature PUBLIC lib${CMAKE_PROJECT_NAME} CLI11::CLI11)
target_link_libraries(clustering-sample PUBLIC lib${CMAKE_PROJECT_NAME} CLI11::CLI11)
target_link_libraries(clustering-secondary-sample PUBLIC lib${CMAKE_PROJECT_NAME} CLI11::CLI11)
target_link_libraries(clustering-solve PUBLIC lib${CMAKE_PROJECT_NAME} CLI11::CLI11)
target_link_libraries(clustering-uncut PUBLIC lib${CMAKE_PROJECT_NAME} CLI11::CLI11)
target_link_libraries(clustering-endend PUBLIC lib${CMAKE_PROJECT_NAME} CLI11::CLI11)
target_link_libraries(strokestrip PUBLIC lib${CMAKE_PROJECT_NAME} CLI11::CLI11)

# Spdlog
include(spdlog)
target_link_libraries(lib${CMAKE_PROJECT_NAME} PUBLIC spdlog::spdlog)

# Serialization
include(nlohmannjson)
target_link_libraries(lib${CMAKE_PROJECT_NAME} PUBLIC nlohmannjson::nlohmannjson)

# Pybind
add_subdirectory(python/core)

# Gurobi
if (GUROBI_DIR STREQUAL "")
    set(GUROBI_DIR "$ENV{GUROBI_HOME}")

    # Ok, I used the 'cmake gui' to build with VS17, which required several cache deletes, 
    # then I used ctrl-shift-p in vscode to do build. 
else()
    message(STATUS "Using GUROBI_DIR: ${GUROBI_DIR}")
endif()

find_package(GUROBI REQUIRED)
include_directories(${GUROBI_INCLUDE_DIRS})
target_include_directories(lib${CMAKE_PROJECT_NAME} PUBLIC ${GUROBI_INCLUDE_DIRS})
target_link_libraries(lib${CMAKE_PROJECT_NAME} PUBLIC ${GUROBI_LIBRARY})
target_link_libraries(lib${CMAKE_PROJECT_NAME} PUBLIC optimized ${GUROBI_CXX_LIBRARY} debug ${GUROBI_CXX_DEBUG_LIBRARY})

# Eigen
# find_package(Eigen3 CONFIG REQUIRED)
include(eigen3)
target_link_libraries(lib${CMAKE_PROJECT_NAME} PUBLIC Eigen3::Eigen)

# TBB
include(tbb)
target_link_libraries(lib${CMAKE_PROJECT_NAME} PUBLIC TBB::tbb)

# span-lite
include(span-lite)
target_link_libraries(lib${CMAKE_PROJECT_NAME} PUBLIC nonstd::span-lite)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    include(FeatureSummary)
    feature_summary(WHAT ALL)
endif()

# Glm
add_subdirectory(external/glm)
target_link_libraries(lib${CMAKE_PROJECT_NAME} PUBLIC glm)

# Cornucopia
add_subdirectory(external/Cornucopia)
target_link_libraries(lib${CMAKE_PROJECT_NAME} PUBLIC Cornucopia)

target_include_directories(lib${CMAKE_PROJECT_NAME} PUBLIC external/)

# SketchConnectivity
add_subdirectory(external/SketchConnectivity)
target_link_libraries(lib${CMAKE_PROJECT_NAME} PUBLIC sketch_connectivity)
