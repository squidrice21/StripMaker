/**
 * @brief This header file contains functions for serializing and deserializing
 * Junction objects to and from JSON format.
 */
#pragma once

#include "../stroke_strip_src/Serialization.h"
#include "Junction.h"

#include <glm/glm.hpp>
#include <nlohmann/json.hpp>

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace {
using json_t = nlohmann::json;
}

void serialize(json_t &j, const Junction &v);
void deserialize(Junction &v, const json_t &j);

void read_junc_json(const std::string &filename,
                    std::vector<Junction> &junctions);
void write_junc_json(const std::string &filename,
                     const std::vector<Junction> &junctions);
