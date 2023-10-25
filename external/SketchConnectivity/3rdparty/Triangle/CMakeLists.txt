project(Triangle)

set(HEADER_FILES Triangle/triangle.h)
set(SOURCE_FILES Triangle/triangle.c)

add_library(${PROJECT_NAME} ${HEADER_FILES} ${SOURCE_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC .)
target_compile_definitions(${PROJECT_NAME} PUBLIC ANSI_DECLARATORS)
target_compile_definitions(${PROJECT_NAME} PRIVATE NO_TIMER TRILIBRARY)
if (UNIX AND NOT APPLE)
    target_compile_options(${PROJECT_NAME} PRIVATE "-fPIC")
endif()
