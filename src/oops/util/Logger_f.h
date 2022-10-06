/*
 * (C) Copyright 2022 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_LOGGER_F_H_
#define OOPS_UTIL_LOGGER_F_H_

#include <cstdint>
#include <cstring>

using int32 = std::int32_t;

namespace oops {

extern "C" {

void log_info_f(char* msg, int32 newl, int32 flush);
void log_error_f(char* msg, int32 newl, int32 flush);
void log_warning_f(char* msg, int32 newl, int32 flush);
void log_debug_f(char* msg, int32 newl, int32 flush);
void log_trace_f(char* msg, int32 newl, int32 flush);
void log_stats_f(char* msg, int32 newl, int32 flush);
void log_test_f(char* msg, int32 newl, int32 flush);

}

}  // namespace oops

#endif  // OOPS_UTIL_LOGGER_F_H_
