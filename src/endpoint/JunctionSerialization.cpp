#include "JunctionSerialization.h"

#include <cstdint>
#include <fstream>

void serialize(json_t &j, const Junction &v) {
  j = json_t::object();
  serialize(j["from_sid"], v.from.first);
  serialize(j["from_time"], v.from.second);
  serialize(j["to_sid"], v.to.first);
  serialize(j["to_time"], v.to.second);
  serialize(j["from_pos"], v.from_pos);
  serialize(j["to_pos"], v.to_pos);
  serialize(j["probability"], v.probability);
  std::uint8_t tmp_type = std::uint8_t(v.type);
  serialize(j["type"], tmp_type);
}
void deserialize(Junction &v, const json_t &j) {
  deserialize(v.from.first, j["from_sid"]);
  deserialize(v.from.second, j["from_time"]);
  deserialize(v.to.first, j["to_sid"]);
  deserialize(v.to.second, j["to_time"]);
  deserialize(v.from_pos, j["from_pos"]);
  deserialize(v.to_pos, j["to_pos"]);
  deserialize(v.probability, j["probability"]);
  std::uint8_t tmp_type;
  deserialize(tmp_type, j["type"]);
  v.type = sketching::JunctionType::Type(tmp_type);
}

void read_junc_json(const std::string &filename,
                    std::vector<Junction> &junctions) {
  std::ifstream fstr(filename);
  json_t json;
  fstr >> json;
  auto &junction_json = json["m_junction"];

  deserialize(junctions, junction_json);
}

void write_junc_json(const std::string &filename,
                     const std::vector<Junction> &junctions) {
  json_t json;

  auto &junction_json = json["m_junction"];
  junction_json = json_t::object();

  serialize(junction_json, junctions);
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
    fstr << json.dump(1, '\t');
  }
}
