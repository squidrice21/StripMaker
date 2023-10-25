if(TARGET TBB::tbb)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    tbb
    GIT_REPOSITORY https://github.com/oneapi-src/oneTBB.git
    GIT_TAG 3df08fe234f23e732a122809b40eb129ae22733f #v2021.5.0
    GIT_SHALLOW FALSE
)

set(TBB_TEST OFF CACHE BOOL "Enable testing" FORCE)

FetchContent_MakeAvailable(tbb)
