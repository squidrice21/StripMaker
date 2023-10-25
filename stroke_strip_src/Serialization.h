#pragma once

#include "Cluster.h"
#include "Fitting.h"

#include <glm/glm.hpp>
#include <nlohmann/json.hpp>

#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace {
using json_t = nlohmann::json;
}

// Serialization of basic types
void serialize(json_t &j, const int &v);
void deserialize(int &v, const json_t &j);
void serialize(json_t &j, const size_t &v);
void deserialize(size_t &v, const json_t &j);
void serialize(json_t &j, const bool &v);
void deserialize(bool &v, const json_t &j);
void serialize(json_t &j, const float &v);
void deserialize(float &v, const json_t &j);
void serialize(json_t &j, const double &v);
void deserialize(double &v, const json_t &j);
void serialize(json_t &j, const std::string &v);
void deserialize(std::string &v, const json_t &j);
void serialize(json_t &j, const std::uint8_t &v);
void deserialize(std::uint8_t &v, const json_t &j);

void serialize(json_t &j, const glm::dvec2 &v);
void deserialize(glm::dvec2 &v, const json_t &j);

void serialize(json_t &j, const Cluster::Stroke &v);
void deserialize(Cluster::Stroke &v, const json_t &j);
void serialize(json_t &j, const Cluster::XSecPoint &v);
void deserialize(Cluster::XSecPoint &v, const json_t &j);
void serialize(json_t &j, const Cluster::XSecConnection &v);
void deserialize(Cluster::XSecConnection &v, const json_t &j);
void serialize(json_t &j, const Cluster::XSec &v);
void deserialize(Cluster::XSec &v, const json_t &j);
void serialize(json_t &j, const FittedCurve &v);
void deserialize(FittedCurve &v, const json_t &j);

// Serialization of array types
template <typename T, typename GetterFunc>
void serialize_array(json_t &v, const GetterFunc &getter,
                     const std::size_t num_values) {
  v = json_t::array();
  for (std::size_t i = 0; i < num_values; ++i) {
    json_t subv;
    T ith_value = getter(i);
    serialize(subv, ith_value);
    v.push_back(subv);
  }
}
template <typename T, typename Allocator>
void serialize(json_t &j, const std::vector<T, Allocator> &vec) {
  serialize_array<T>(
    j, [&vec](const std::size_t i) { return vec[i]; }, vec.size());
}

template <typename T, typename ResizeFunc, typename SetValueFunc>
void deserialize_array(const json_t &v, const ResizeFunc &resize,
                       const SetValueFunc &set_value) {
  resize(v.size());
  for (std::size_t i = 0; i < v.size(); ++i) {
    T v_;
    deserialize(v_, v[i]);
    set_value(i, v_);
  }
}
template <typename T, typename Allocator>
void deserialize(std::vector<T, Allocator> &vec, const json_t &v) {
  deserialize_array<T>(
    v, [&vec](const std::size_t sz) { vec.resize(sz); },
    [&vec](const std::size_t i, const T &v_) { vec[i] = v_; });
}

void read_json(const std::string &filename, Cluster &cluster);
void write_json(const std::string &filename, const Cluster &cluster);
