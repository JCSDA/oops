/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_STATEINFO_H_
#define OOPS_BASE_STATEINFO_H_

#include <string>

#include "oops/base/PostBase.h"
#include "oops/util/Logger.h"

namespace eckit {
  class Configuration;
}

namespace oops {

/// Handles writing-out of forecast fields.
/*!
 *  Write out forecast fields.
 */

template <typename FLDS> class StateInfo : public PostBase<FLDS> {
 public:
  StateInfo(const std::string sgrep, const PostTimerParameters & timerParameters):
    PostBase<FLDS>(timerParameters), sgrep_(sgrep) {}
  StateInfo(const std::string sgrep, const eckit::Configuration & conf):
    StateInfo(sgrep, validateAndDeserialize<PostTimerParameters>(conf)) {}
  StateInfo(const std::string sgrep, const util::Duration & freq):
    PostBase<FLDS>(freq), sgrep_(sgrep) {}
  ~StateInfo() {}

 private:
  const std::string sgrep_;
  void doProcessing(const FLDS & xx) override {Log::info() << sgrep_ << xx << std::endl;}
};

}  // namespace oops

#endif  // OOPS_BASE_STATEINFO_H_
