#include "Logger.h"
#include "spdlog/common.h"

#include <spdlog/sinks/stdout_color_sinks.h>

// Retrieve current logger
spdlog::logger &logger() {
  static auto default_logger = spdlog::stdout_color_mt("Clustering");
  default_logger->set_level(spdlog::level::warn);
  return *default_logger;
}
