project(sketching-python)

# pybind11
include(pybind11)
pybind11_add_module(${PROJECT_NAME} ${HEADER_FILES} ${SOURCE_FILES})

file(GLOB_RECURSE HEADER_FILES CONFIGURE_DEPENDS *.h)
file(GLOB_RECURSE SOURCE_FILES CONFIGURE_DEPENDS *.cpp)
list(FILTER SOURCE_FILES EXCLUDE REGEX ".*Main.*.cpp$")
list(FILTER HEADER_FILES EXCLUDE REGEX ".*Main.*.h$")

# So we can do `import _sketching` from Python.
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME _sketching)
target_link_libraries(${PROJECT_NAME} PRIVATE libstripmaker)
target_link_libraries(${PROJECT_NAME} PUBLIC pybind11::pybind11)
