if(TARGET spdlog::spdlog)
  return()
endif()

include(FetchContent)
  FetchContent_Declare(
  spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG        v1.5.0
)
message(STATUS "creating target 'spdlog::spdlog'")

FetchContent_MakeAvailable(spdlog)

set_target_properties(spdlog PROPERTIES FOLDER third_party)

## Kill extra warnings
if(
  "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" OR
  "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"
)
target_compile_options(spdlog PRIVATE
    "-Wno-sign-conversion"
)
endif()