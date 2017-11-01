/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef UTIL_LOGGER_H_
#define UTIL_LOGGER_H_

#include "eckit/log/Log.h"
#include "util/LibOOPS.h"

namespace oops {

// -----------------------------------------------------------------------------

struct Log {
    static std::ostream& info()    { return eckit::Log::info(); }
    static std::ostream& error()   { return eckit::Log::error(); }
    static std::ostream& warning() { return eckit::Log::warning(); }
    static std::ostream& debug()   { return eckit::Log::debug<LibOOPS>(); }

    // Following are non-default to eckit. They wrap eckit::Log::info() with additional prefix
    static std::ostream& trace()   { return LibOOPS::instance().traceChannel(); }  // prefix "OOPS_TRACE"
    static std::ostream& stats()   { return LibOOPS::instance().statsChannel(); }  // prefix "OOPS_STATS"
    static std::ostream& test()    { return LibOOPS::instance().testChannel(); }   // prefix "Test     :"
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // UTIL_LOGGER_H_
