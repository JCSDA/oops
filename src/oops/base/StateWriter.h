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
#include "oops/base/PostTimerParameters.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/parameters/Parameters.h"

namespace oops {

template <typename FLDS> class StateWriterParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(StateWriterParameters, Parameters)

 public:
  /// \brief Options determining the time steps at which the state is written out.
  PostTimerParameters postTimer{this};
  /// \brief Options passed to the FLDS::write() function.
  typename FLDS::WriteParameters_ write{this};
};

/// Handles writing-out of forecast fields.
/*!
 *  Write out forecast fields.
 */

template <typename FLDS> class StateWriter : public PostBase<FLDS> {
 public:
  explicit StateWriter(const StateWriterParameters<FLDS> & parameters):
    PostBase<FLDS>(parameters.postTimer),
    writeParameters_(parameters.write) {}
  explicit StateWriter(const eckit::Configuration & conf):
    // NOLINTNEXTLINE(runtime/explicit): lint misinterprets the next line as an implicit constructor
    StateWriter(validateAndDeserialize<StateWriterParameters<FLDS>>(conf)) {}
  ~StateWriter() {}

 private:
  const typename FLDS::WriteParameters_ writeParameters_;
  void doProcessing(const FLDS & xx) override {xx.write(writeParameters_);}
};

}  // namespace oops

#endif  // OOPS_BASE_STATEWRITER_H_
