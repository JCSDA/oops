/*
 * (C) Copyright 2022 MetOffice.
 * (C) Copyright 2021-2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_HYBRIDLINEARMODELCOEFFS_H_
#define OOPS_GENERIC_HYBRIDLINEARMODELCOEFFS_H_

#include <map>
#include <memory>
#include <string>

#include "eckit/mpi/Comm.h"
#include "oops/base/Geometry.h"
#include "oops/base/Variables.h"
#include "oops/generic/HtlmEnsemble.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {
template <typename MODEL>
class HybridLinearModelCoeffsParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(HybridLinearModelCoeffsParameters, Parameters)
  typedef HtlmEnsembleParameters<MODEL>      HtlmEnsembleParameters_;
  typedef typename Geometry<MODEL>::Parameters_ GeometryParameters_;

 public:
  RequiredParameter< HtlmEnsembleParameters_> htlmEnsemble{"htlm ensemble", this};
  /// Forecast length.
  RequiredParameter<util::Duration> windowLength{"window length", this};
  /// Window begin
  RequiredParameter<util::DateTime> windowBegin{"window begin", this};
};

template <typename MODEL>
class HybridLinearModelCoeffs{
 public:
  typedef HybridLinearModelCoeffsParameters<MODEL>      HybridLinearModelCoeffsParameters_;
  typedef HtlmEnsemble<MODEL>                           HtlmEnsemble_;
  typedef HtlmEnsembleParameters<MODEL>                 HtlmEnsembleParameters_;
  typedef Geometry<MODEL>                               Geometry_;

  static const std::string classname() {return "oops::HybridLinearCoeffs";}

  /// constructor
  /* The geometry for the state is a yaml parameter
   The geometry passed in is used for the TLM geometry in HtlmEnsemble
   This is ultimatly first constucted in HybridLinearModel fromits paramters */
  HybridLinearModelCoeffs(const HybridLinearModelCoeffsParameters_ &, const Geometry_ &,
                                                                  const util::Duration &);

 private:
  const HybridLinearModelCoeffsParameters_  params_;
  HtlmEnsemble_ ens_;
  const util::DateTime windowBegin_;
  const util::Duration windowLength_;
  const util::Duration tstep_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
HybridLinearModelCoeffs<MODEL>::HybridLinearModelCoeffs(const  HybridLinearModelCoeffsParameters_
                                                        & params,  const Geometry_ & geomTLM,
                                                        const util::Duration & tstep)
: params_(params),  ens_(params.htlmEnsemble.value(), geomTLM),
  windowBegin_(params_.windowBegin.value()),
  windowLength_(params_.windowLength.value()), tstep_(tstep) {
  util::DateTime time(windowBegin_);
  /// step ensemble
  while (time < (windowBegin_ + windowLength_)) {
    time += tstep_;
    ens_.step(tstep_);
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_HYBRIDLINEARMODELCOEFFS_H_
