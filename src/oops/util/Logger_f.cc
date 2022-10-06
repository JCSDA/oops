/*
 * (C) Copyright 2022 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/util/Logger.h"
#include "oops/util/Logger_f.h"

using int32 = std::int32_t;

namespace oops {

void log_generic(std::ostream * stream, char* msg, int32 newl, int32 flush) {
    if (::strlen(msg)) {
        *stream << msg;
    } else {
        *stream << " ";
    }
    if (newl) *stream << eckit::newl;
    if (flush) *stream << std::flush;
}
// -----------------------------------------------------------------------------
void log_info_f(char* msg, int32 newl, int32 flush) {
    log_generic(&Log::info(), msg, newl, flush);
}
// -----------------------------------------------------------------------------
void log_error_f(char* msg, int32 newl, int32 flush) {
    log_generic(&Log::error(), msg, newl, flush);
}
// -----------------------------------------------------------------------------
void log_warning_f(char* msg, int32 newl, int32 flush) {
    log_generic(&Log::warning(), msg, newl, flush);
}
// -----------------------------------------------------------------------------
void log_debug_f(char* msg, int32 newl, int32 flush) {
    log_generic(&Log::debug(), msg, newl, flush);
}
// -----------------------------------------------------------------------------
void log_trace_f(char* msg, int32 newl, int32 flush) {
    log_generic(&Log::trace(), msg, newl, flush);
}
// -----------------------------------------------------------------------------
void log_stats_f(char* msg, int32 newl, int32 flush) {
    log_generic(&Log::stats(), msg, newl, flush);
}
// -----------------------------------------------------------------------------
void log_test_f(char* msg, int32 newl, int32 flush) {
    log_generic(&Log::test(), msg, newl, flush);
}

}  // namespace oops
