#include "Serialization.h"

#include "Cluster.h"

#include <fstream>
#include <regex>

void serialize(json_t &j, const int &v) { j = v; }
void deserialize(int &v, const json_t &j) { v = j; }

void serialize(json_t &j, const size_t &v) { j = v; }
void deserialize(size_t &v, const json_t &j) { v = j; }

void serialize(json_t &j, const bool &v) { j = v; }
void deserialize(bool &v, const json_t &j) { v = j; }

void serialize(json_t &j, const float &v) { j = v; }
void deserialize(float &v, const json_t &j) { v = j; }

void serialize(json_t &j, const double &v) { j = v; }
void deserialize(double &v, const json_t &j) { v = j; }

void serialize(json_t &j, const std::string &v) { j = v; }
void deserialize(std::string &v, const json_t &j) { v = j; }

void serialize(json_t &j, const std::uint8_t &v) { j = v; }
void deserialize(std::uint8_t &v, const json_t &j) { v = j; }

void serialize(json_t &j, const glm::dvec2 &v) {
  j = json_t::array();
  json_t v1, v2;
  serialize(v1, v.x);
  serialize(v2, v.y);
  j.push_back(v1);
  j.push_back(v2);
}
void deserialize(glm::dvec2 &v, const json_t &j) {
  deserialize(v.x, j[0]);
  deserialize(v.y, j[1]);
}

void serialize(json_t &j, const Cluster::Stroke &v) {
  j = json_t::object();
  serialize(j["stroke_ind"], v.stroke_ind);
  serialize(j["cluster_ind"], v.cluster_ind);
  serialize<glm::dvec2>(j["points"], v.points);
  serialize<double>(j["u"], v.u);
  serialize(j["spiral"], v.spiral);
}
void deserialize(Cluster::Stroke &v, const json_t &j) {
  deserialize(v.stroke_ind, j["stroke_ind"]);
  deserialize(v.cluster_ind, j["cluster_ind"]);
  deserialize<glm::dvec2>(v.points, j["points"]);
  deserialize<double>(v.u, j["u"]);
  if (j.contains("spiral"))
    deserialize(v.spiral, j["spiral"]);
}

void serialize(json_t &j, const Cluster::XSecPoint &v) {
  serialize(j["stroke_ind"], v.stroke_ind);
  serialize(j["cluster_ind"], v.cluster_ind);
  serialize(j["stroke_idx"], v.stroke_idx);
  serialize(j["i"], v.i);
  serialize(j["point"], v.point);
  serialize(j["to_next"], v.to_next);
  serialize(j["tangent"], v.tangent);
  serialize(j["stroke_within_idx"], v.stroke_within_idx);
  serialize(j["stroke_xsec_count"], v.stroke_xsec_count);
}
void deserialize(Cluster::XSecPoint &v, const json_t &j) {
  deserialize(v.stroke_ind, j["stroke_ind"]);
  deserialize(v.cluster_ind, j["cluster_ind"]);
  deserialize(v.stroke_idx, j["stroke_idx"]);
  deserialize(v.i, j["i"]);
  deserialize(v.point, j["point"]);
  deserialize(v.to_next, j["to_next"]);
  deserialize(v.tangent, j["tangent"]);
  if (j.contains("stroke_within_idx"))
    deserialize(v.stroke_within_idx, j["stroke_within_idx"]);
  if (j.contains("stroke_xsec_count"))
    deserialize(v.stroke_xsec_count, j["stroke_xsec_count"]);
}

void serialize(json_t &j, const Cluster::XSecConnection &v) {
  serialize(j["a_idx"], v.a_idx);
  serialize(j["b_idx"], v.b_idx);
  serialize(j["weight"], v.weight);
}
void deserialize(Cluster::XSecConnection &v, const json_t &j) {
  deserialize(v.a_idx, j["a_idx"]);
  deserialize(v.b_idx, j["b_idx"]);
  deserialize(v.weight, j["weight"]);
}

void serialize(json_t &j, const Cluster::XSec &v) {
  serialize<Cluster::XSecPoint>(j["points"], v.points);
  serialize<Cluster::XSecConnection>(j["connections"], v.connections);
  serialize(j["center_idx"], v.center_idx);
  serialize(j["u"], v.u);
  serialize(j["connector"], v.connector);
}
void deserialize(Cluster::XSec &v, const json_t &j) {
  deserialize<Cluster::XSecPoint>(v.points, j["points"]);
  deserialize<Cluster::XSecConnection>(v.connections, j["connections"]);
  deserialize(v.center_idx, j["center_idx"]);
  deserialize(v.u, j["u"]);
  deserialize(v.connector, j["connector"]);
}

void serialize(json_t &j, const FittedCurve &v) {
  serialize<glm::dvec2>(j["centerline"], v.centerline);
  serialize<double>(j["widths"], v.widths);
  serialize(j["fit_sample_xsecs"], v.fit_sample_xsecs);
  serialize(j["cluster_idx"], v.cluster_idx);
}
void deserialize(FittedCurve &v, const json_t &j) {
  deserialize<glm::dvec2>(v.centerline, j["centerline"]);
  deserialize<double>(v.widths, j["widths"]);
  if (j.contains("fit_sample_xsecs"))
    deserialize(v.fit_sample_xsecs, j["fit_sample_xsecs"]);
  if (j.contains("cluster_idx"))
    deserialize(v.cluster_idx, j["cluster_idx"]);
}

void read_json(const std::string &filename, Cluster &cluster) {
  std::ifstream fstr(filename);
  json_t json;
  fstr >> json;
  auto &cluster_json = json["m_cluster"];

  deserialize(cluster.strokes, cluster_json["properties"]["strokes"]);
  if (cluster_json["properties"].contains("original_input_strokes"))
    deserialize(cluster.original_input_strokes,
                cluster_json["properties"]["original_input_strokes"]);
  deserialize(cluster.xsecs, cluster_json["properties"]["xsecs"]);
  deserialize(cluster.periodic, cluster_json["properties"]["periodic"]);
  deserialize(cluster.fit, cluster_json["properties"]["fit"]);

  for (auto &k_v : cluster_json["properties"].items()) {
    std::string k = k_v.key();
    if (k.find("obj_") != std::string::npos) {
      std::string obj_k = std::regex_replace(k, std::regex("obj_"), "");
      double v;
      deserialize(v, k_v.value());
      cluster.obj_term_values[obj_k] = v;
    }
  }

  // Verify the file
  std::string stroke_str;
  for (auto const &s : cluster.strokes) {
    if (!stroke_str.empty())
      stroke_str += ' ';
    stroke_str += std::to_string(s.stroke_ind);
  }
  std::string in_stroke_str;
  deserialize(in_stroke_str, json["clusters_str"]);
  assert(in_stroke_str == stroke_str);
}

void write_json(const std::string &filename, const Cluster &cluster) {
  json_t json;

  std::string stroke_str;
  for (auto const &s : cluster.strokes) {
    if (!stroke_str.empty())
      stroke_str += ' ';
    stroke_str += std::to_string(s.stroke_ind);
  }
  serialize(json["clusters_str"], stroke_str);

  auto &cluster_json = json["m_cluster"];
  cluster_json = json_t::object();
  cluster_json["properties"] = json_t::object();

  serialize(cluster_json["properties"]["strokes"], cluster.strokes);
  serialize(cluster_json["properties"]["original_input_strokes"],
            cluster.original_input_strokes);
  serialize(cluster_json["properties"]["xsecs"], cluster.xsecs);
  // serialize(cluster_json["properties"]["potential_strokes"],
  // cluster.potential_strokes);
  serialize(cluster_json["properties"]["periodic"], cluster.periodic);
  serialize(cluster_json["properties"]["fit"], cluster.fit);

  // Parameterization related terms
  for (const auto &term_value : cluster.obj_term_values)
    serialize(cluster_json["properties"]["obj_" + term_value.first],
              term_value.second);

  // Timestamp
  {
    time_t rawtime;
    struct tm *timeinfo;
    char buffer[80];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer, sizeof(buffer), "%d-%m-%Y %H:%M:%S", timeinfo);
    json["saved_at"] = buffer;
  }

  // Write to file
  {
    std::ofstream fstr(filename);
    // std::cout << "Writing json file: " << filename << std::endl;
    fstr << json.dump(1, '\t');
  }
}
