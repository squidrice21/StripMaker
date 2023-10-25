if(TARGET pybind11::pybind11)
  return()
endif()

include(FetchContent)
  FetchContent_Declare(
    pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.6.2
)
message(STATUS "creating target 'pybind11::pybind11'")

set(BUILD_TESTING OFF CACHE BOOL "Enable the tests" FORCE)

FetchContent_MakeAvailable(pybind11)
