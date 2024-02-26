/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_OBSSPACEBASE_H_
#define OOPS_BASE_OBSSPACEBASE_H_

#include <boost/noncopyable.hpp>

#include "eckit/mpi/Comm.h"

#include "oops/util/DateTime.h"
#include "oops/util/Printable.h"
#include "oops/util/TimeWindow.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// Base class for observation spaces.

class ObsSpaceBase : public util::Printable,
                     private boost::noncopyable {
 public:
  ObsSpaceBase(const eckit::Configuration &, const eckit::mpi::Comm &,
               const util::TimeWindow &);

  virtual ~ObsSpaceBase() {}

/// Access information
  util::DateTime windowStart() const {return timeWindow_.start();}
  util::DateTime windowEnd() const {return timeWindow_.end();}
  util::TimeWindow timeWindow() const {return timeWindow_;}

  int64_t getSeed() const {return seed_;}

 private:
  static int instances_;

  const util::TimeWindow timeWindow_;

  int64_t seed_;
  const int instance_;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSSPACEBASE_H_
