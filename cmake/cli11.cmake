if(TARGET CLI11::CLI11)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    cli11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG v2.1.2
)
message(STATUS "creating target 'CLI11::CLI11'")
FetchContent_MakeAvailable(cli11)
