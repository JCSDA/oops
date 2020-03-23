/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_STATEWRITER_H_
#define OOPS_BASE_STATEWRITER_H_

#include "oops/base/PostBase.h"

#include "eckit/config/LocalConfiguration.h"
#include "oops/util/DateTime.h"

namespace oops {

/// Handles writing-out of forecast fields.
/*!
 *  Write out forecast fields.
 */

template <typename FLDS> class StateWriter : public PostBase<FLDS> {
 public:
  explicit StateWriter(const eckit::Configuration & conf):
    PostBase<FLDS>(conf), ppConfig_(conf) {}
  ~StateWriter() {}

 private:
  const eckit::LocalConfiguration ppConfig_;
  void doProcessing(const FLDS & xx) override {xx.write(ppConfig_);}
};

}  // namespace oops

#endif  // OOPS_BASE_STATEWRITER_H_
