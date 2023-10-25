#pragma once

#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>

/**
Retrieve the logger

@return     A const reference to logger object.

Usage:

- Logging something

  There are the following functions that allow logging with different severities.
  If the severity of the message is bigger than the logger's, the message is printed.
  Here are the severities in order of low importance to high.

  trace, debug, info, warn, err, critical

  Example:
  logger().trace("I have {} apples", 3);

- Setting the logger() severity:

  Default severity is info, you can change it like this:

  logger().set_severity(spdlog::level::trace);

- Adding a sink (where logs are printed):

  By default only one sink is present that prints to the screen.
  By adding this to the top of your main file you can also start printing
  to a file.

  #include <spdlog/sinks/basic_file_sink.h>
  logger().sinks().push_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>("filename.txt"));
*/
spdlog::logger& logger();
