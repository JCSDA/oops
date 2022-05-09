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

#include <memory>
#include <sstream>
#include <string>

#include "eckit/system/Library.h"
#include "eckit/utils/Translator.h"

#include "oops/util/TestReference.h"

namespace oops {

// -----------------------------------------------------------------------------

class LibOOPS : public eckit::system::Library {
 public:
  LibOOPS();

  ~LibOOPS();

  static LibOOPS& instance();

  eckit::Channel& infoChannel() const;
  eckit::Channel& debugChannel() const override;

  eckit::Channel& traceChannel() const;
  eckit::Channel& statsChannel() const;
  eckit::Channel& testChannel() const;
  eckit::Channel& timerChannel() const;

  void initialise();
  void testReferenceInitialise(const eckit::LocalConfiguration &);
  void teeOutput(const std::string &);
  void finalise(bool finaliseMPI = true);

 protected:
  const void* addr() const override;

  std::string version() const override;

  std::string gitsha1(unsigned int count) const override;

  mutable std::unique_ptr<eckit::Channel> infoChannel_;
  mutable std::unique_ptr<eckit::Channel> debugChannel_;

  mutable std::unique_ptr<eckit::Channel> traceChannel_;
  mutable std::unique_ptr<eckit::Channel> statsChannel_;
  mutable std::unique_ptr<eckit::Channel> testChannel_;
  mutable std::unique_ptr<eckit::Channel> timerChannel_;

  size_t rank_;
  bool debug_;
  std::string predebug_;
  bool trace_;
  std::string pretrace_;

  // TestReferece associated member variables
  std::stringstream testStream_;
  TestReference testReference_;

 private:
  bool enable_timer_channel_;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_UTIL_LIBOOPS_H_
