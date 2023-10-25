project(GraphicsGems)

set(HEADER_FILES
    src/GraphicsGems/GraphicsGems.h
    src/GraphicsGems/FitCurves.h
)
set(SOURCE_FILES
    src/GraphicsGems/GraphicsGems.c
    src/GraphicsGems/FitCurves.c
)

add_library(${PROJECT_NAME} ${HEADER_FILES} ${SOURCE_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC src)
if (UNIX AND NOT APPLE)
    target_compile_options(${PROJECT_NAME} PRIVATE "-fPIC")
endif()
