cmake_minimum_required(VERSION 3.15)

project(sketch_connectivity)

option(WITH_GUROBI "Include features which depend on Gurobi (requires Gurobi to be installed)" ON)
option(SPDLOG_HEADER_ONLY "Use the header-only version of spdlog" OFF)
option(ASAN "Enable address sanitizer" OFF)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

file(GLOB_RECURSE HEADER_FILES CONFIGURE_DEPENDS *.h)
file(GLOB_RECURSE SOURCE_FILES CONFIGURE_DEPENDS *.cpp)

if (MSVC)
    string(REPLACE "/Ob0" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
endif()

add_library(${PROJECT_NAME}
            ${HEADER_FILES} ${SOURCE_FILES} ${NATVIS_FILES})
source_group(TREE ${PROJECT_SOURCE_DIR}/ FILES ${HEADER_FILES} ${SOURCE_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC ./)
set_property(TARGET ${PROJECT_NAME} PROPERTY PROJECT_LABEL lib${PROJECT_NAME})
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)
target_compile_definitions(${PROJECT_NAME} PRIVATE _USE_MATH_DEFINES)
if (MSVC)
    target_compile_options(
        ${PROJECT_NAME} PRIVATE
        /W4 # Warning level
        /we4834 # ignoring nodiscard is an error
        )
    target_compile_definitions(${PROJECT_NAME} PRIVATE
                               _CRT_SECURE_NO_WARNINGS
                               _SILENCE_CXX17_OLD_ALLOCATOR_MEMBERS_DEPRECATION_WARNING)
    target_compile_options(${PROJECT_NAME} PRIVATE $<$<CONFIG:DEBUG>:/Ob1>)
    if (ASAN)
        target_compile_options(${PROJECT_NAME} PRIVATE /fsanitize=address)
    endif()
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -pedantic -Wno-ignored-qualifiers)
    if (ASAN)
        target_compile_options(${PROJECT_NAME} PRIVATE -fsanitize=address)
    endif()
endif()
if (UNIX AND NOT APPLE)
    target_compile_options(${PROJECT_NAME} PRIVATE -fPIC)
endif()
target_precompile_headers(${PROJECT_NAME} PRIVATE
    <Eigen/Core>
    <nonstd/span.hpp>
    <spdlog/spdlog.h>
    <algorithm>
    <map>
    <memory>
    <set>
    <string>
    <unordered_map>
    <vector>)

find_package(Eigen3 3.3.7 CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)

# find_package(fmt CONFIG REQUIRED)
# target_link_libraries(${PROJECT_NAME} PUBLIC fmt::fmt)

find_package(spdlog CONFIG REQUIRED)
if (SPDLOG_HEADER_ONLY)
    target_link_libraries(${PROJECT_NAME} PUBLIC spdlog::spdlog_header_only)
else()
    target_link_libraries(${PROJECT_NAME} PUBLIC spdlog::spdlog)
endif()

add_subdirectory(3rdparty/GraphicsGems)
set_target_properties(GraphicsGems PROPERTIES FOLDER "3rdparty")
target_link_libraries(${PROJECT_NAME} PUBLIC GraphicsGems)

add_subdirectory(3rdparty/span-lite)
target_link_libraries(${PROJECT_NAME} PUBLIC span_lite)
target_compile_definitions(${PROJECT_NAME} PUBLIC span_CONFIG_NO_EXCEPTIONS)

add_subdirectory(3rdparty/Triangle)
set_target_properties(Triangle PROPERTIES FOLDER "3rdparty")
target_link_libraries(${PROJECT_NAME} PUBLIC Triangle)
