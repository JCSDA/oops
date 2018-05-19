/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Baudouin Raoult
/// @author Tiago Quintino
/// @date   December 2016

#ifndef OOPS_UTIL_LIBOOPS_H_
#define OOPS_UTIL_LIBOOPS_H_

#include <string>
#include "eckit/system/Library.h"

namespace oops {

// -----------------------------------------------------------------------------

class LibOOPS : public eckit::system::Library {
 public:
  LibOOPS();

  ~LibOOPS();

  static LibOOPS& instance();

  virtual eckit::Channel& traceChannel() const;

  virtual eckit::Channel& statsChannel() const;

  virtual eckit::Channel& testChannel() const;

  void finalise();

 protected:
  const void* addr() const;

  virtual std::string version() const;

  virtual std::string gitsha1(unsigned int count) const;

  mutable eckit::ScopedPtr<eckit::Channel> traceChannel_;

  mutable eckit::ScopedPtr<eckit::Channel> statsChannel_;

  mutable eckit::ScopedPtr<eckit::Channel> testChannel_;

  bool trace_;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_UTIL_LIBOOPS_H_
