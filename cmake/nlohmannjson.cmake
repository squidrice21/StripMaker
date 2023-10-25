
if(TARGET nlohmannjson::nlohmannjson)
  return()
endif()

message(STATUS "Creating target 'nlohmannjson::nlohmannjson'")

# nlohmannjson is a big repo for a single header, so we just download the release archive
set(NLOHMANNJSON_VERSION "v3.7.3")

include(FetchContent)
FetchContent_Declare(
  nlohmannjson
  URL "https://github.com/nlohmann/json/releases/download/${NLOHMANNJSON_VERSION}/include.zip"
  URL_HASH SHA256=87b5884741427220d3a33df1363ae0e8b898099fbc59f1c451113f6732891014
)
FetchContent_MakeAvailable(nlohmannjson)

add_library(nlohmannjson::nlohmannjson INTERFACE IMPORTED GLOBAL)
target_include_directories(nlohmannjson::nlohmannjson INTERFACE ${nlohmannjson_SOURCE_DIR}/include)
