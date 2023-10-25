if(TARGET nonstd::span-lite)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    span_lite
    GIT_REPOSITORY https://github.com/martinmoene/span-lite.git
    GIT_TAG v0.10.3
)
message(STATUS "creating target 'nonstd::span-lite'")
FetchContent_MakeAvailable(span_lite)
