if(TARGET Eigen3::Eigen)
  return()
endif()

include(FetchContent)
  FetchContent_Declare(
    eigen3
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG        3.4.0
)
message(STATUS "creating target 'Eigen3::Eigen'")

FetchContent_MakeAvailable(eigen3)
