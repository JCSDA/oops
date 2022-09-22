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
#include <vector>

#include "eckit/mpi/Comm.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/Variables.h"
#include "oops/generic/HtlmCalculator.h"
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
  typedef HtlmCalculatorParameters<MODEL>    HtlmCalculatorParameters_;
  typedef typename Geometry<MODEL>::Parameters_ GeometryParameters_;

 public:
  // Ensemble parameters
  RequiredParameter< HtlmEnsembleParameters_> htlmEnsemble{"htlm ensemble", this};
  // Coefficent calculator parameters
  RequiredParameter<HtlmCalculatorParameters_> htlmCalculator{"htlm calculator", this};
  // Forecast length.
  RequiredParameter<util::Duration> windowLength{"window length", this};
  // Window begin
  RequiredParameter<util::DateTime> windowBegin{"window begin", this};
  // variables
  RequiredParameter<Variables> vars{"variables", this};
};

template <typename MODEL>
class HybridLinearModelCoeffs{
 public:
  typedef HybridLinearModelCoeffsParameters<MODEL>      HybridLinearModelCoeffsParameters_;
  typedef HtlmEnsemble<MODEL>                           HtlmEnsemble_;
  typedef HtlmEnsembleParameters<MODEL>                 HtlmEnsembleParameters_;
  typedef HtlmCalculator<MODEL>                         HtlmCalculator_;
  typedef HtlmCalculatorParameters<MODEL>               HtlmCalculatorParameters_;
  typedef Geometry<MODEL>                               Geometry_;
  typedef Increment<MODEL>                              Increment_;

  static const std::string classname() {return "oops::HybridLinearCoeffs";}

  // constructor
  /* The geometry for the state is a yaml parameter
   The geometry passed in is used for the TLM geometry in HtlmEnsemble
   This is ultimately first constructed in HybridLinearModel from its parameters */
  HybridLinearModelCoeffs(const HybridLinearModelCoeffsParameters_ &, const Geometry_ &,
                                                                  const util::Duration &);

 public:
  void updateIncTL(Increment_ &) const;
  void updateIncAD(Increment_ &) const;

 private:
  std::vector<double> getInfluenceVec(const atlas::FieldSet &,
                                        const atlas::idx_t & i,
                                        const atlas::idx_t & k) const;

 private:
  std::map<util::DateTime, atlas::FieldSet> coeffSaver_;
  const Variables vars_;
  const size_t influenceSize_;
  const size_t halfInfluenceSize_;
};

//------------------------------------------------------------------------------

template<typename MODEL>
HybridLinearModelCoeffs<MODEL>::HybridLinearModelCoeffs(const  HybridLinearModelCoeffsParameters_
                                                        & params,  const Geometry_ & geomTLM,
                                                        const util::Duration & tstep) :
    vars_(params.vars.value()), influenceSize_(params.htlmCalculator.value().influenceRegionSize),
    halfInfluenceSize_(influenceSize_/2) {
  HtlmEnsemble_ ens(params.htlmEnsemble.value(), geomTLM);
  HtlmCalculator_ calculator(params.htlmCalculator.value(), params.vars,
    params.htlmEnsemble.value().ensembleSize.value(), geomTLM, params.windowBegin.value());
  const util::DateTime  windowBegin = params.windowBegin.value();
  const util::Duration windowLength = params.windowLength.value();

  util::DateTime time(windowBegin);
  const size_t ensembleSize = params.htlmEnsemble.value().ensembleSize.value();
  // step ensemble and calculate coefficients at time t
  while (time < (windowBegin + windowLength)) {
    time += tstep;
    // step ensemble HtlmEnsemble.h
    ens.step(tstep);
    // create an empty field set to store coefficients at this time
    // dimensions of the set are defined in the calculator object
    // to which it is passed by reference
    atlas::FieldSet coeffFieldSet;
    // emplace coeff field set in map for storage
    coeffSaver_.emplace(time, coeffFieldSet);
    // send ensemble info for coeff calculation and FeildSet for storage of coefficients
    // HtlmCalculator.h
    calculator.calcCoeffs(ens.getLinearEns(), ens.getLinearErrDe(), coeffFieldSet);
    // update increments with coefficients before next step
    for (size_t ensInd = 0; ensInd < ensembleSize; ++ensInd) {
      updateIncTL(ens.getLinearEns()[ensInd]);
    }
  }
}

//------------------------------------------------------------------------------
// updates grid point dx to dx' via dx' = dx + dot(coeffvec,dxinfluenceregion)
template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::updateIncTL(Increment_ & dx) const {
  Log::trace() << "HybridLinearModelCoeffs<MODEL::updateIncTL() starting" << std::endl;
  atlas::FieldSet & dxFset = dx.fieldSet();
  for (size_t varInd = 0; varInd < vars_.size(); ++varInd) {
    auto dxView = atlas::array::make_view<double, 2>(dxFset[varInd]);
    for (atlas::idx_t i = 0; i < dxView.shape(1); ++i) {
      for (atlas::idx_t k = 0; k < dxView.shape(2); ++k) {
        std::vector<double> dxVec = getInfluenceVec(dxFset, i, k);
        auto coeffView = atlas::array::make_view<double, 3>
                                        (coeffSaver_.at(dx.validTime())[vars_[varInd]]);
        for (atlas::idx_t infInd = 0; infInd < atlas::idx_t(dxVec.size()); ++infInd) {
          dxView(i, k) += coeffView(i, k, infInd)*dxVec[infInd];
        }
      }
    }
  }
  Log::trace() << "HybridLinearModelCoeffs<MODEL::updateIncTL() done" << std::endl;
}

//------------------------------------------------------------------------------
// updates grid point dx' to dx via dx = dx - dot(coeffvec,dxinfluenceregion) (adjoint)
template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::updateIncAD(Increment_ & dx) const {
  Log::trace() << "HybridLinearModelCoeffs<MODEL::updateIncAD() starting" << std::endl;
    atlas::FieldSet & dxFset = dx.fieldSet();
  for (size_t varInd = 0; varInd < vars_.size(); ++varInd) {
    auto dxView = atlas::array::make_view<double, 2>(dxFset[varInd]);
    for (atlas::idx_t i = 0; i < dxView.shape(1); ++i) {
      for (atlas::idx_t k = 0; k < dxView.shape(2); ++k) {
        std::vector<double> dxVec = getInfluenceVec(dxFset, i, k);
        auto coeffView = atlas::array::make_view<double, 3>
                                        (coeffSaver_.at(dx.validTime())[vars_[varInd]]);
        for (atlas::idx_t infInd = 0; infInd < atlas::idx_t(dxVec.size()); ++infInd) {
          dxView(i, k) -= coeffView(i, k, infInd)*dxVec[infInd];
        }
      }
    }
  }
  Log::trace() << "HybridLinearModelCoeffs<MODEL::updateIncAD() done" << std::endl;
}

//------------------------------------------------------------------------------
// creates influence region vector for grid point dx, gets points above and below for all vars
template<typename MODEL>
std::vector<double> HybridLinearModelCoeffs<MODEL>::getInfluenceVec(const atlas::FieldSet & dxFset,
                                                           const atlas::idx_t & i,
                                                           const atlas::idx_t & k) const {
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::getInfluenceVec() starting" << std::endl;
  std::vector<double> dxVec(influenceSize_*vars_.size());
  for (size_t varInd = 0; varInd < vars_.size(); ++varInd) {
    const auto dxFview = atlas::array::make_view<double, 2>(dxFset[varInd]);
      for (atlas::idx_t infInd = 0; infInd < influenceSize_; ++infInd) {
        if (k-halfInfluenceSize_ > 0 && k + halfInfluenceSize_ < dxFview.shape(1)) {
          // middle case
          dxVec[vars_.size()*varInd + infInd] = dxFview(i, k-halfInfluenceSize_ + infInd);
        } else if (k-halfInfluenceSize_ <= 0) {
          // start of increment edge case
          dxVec[vars_.size()*varInd + infInd] = dxFview(i, infInd);
        } else if (k+halfInfluenceSize_ >= dxFview.shape(1)) {
          // end of increment edge case
          dxVec[vars_.size()*varInd + infInd] =
                              dxFview(i, (dxFview.shape(1)-influenceSize_) + infInd);
        } else {
          // Checking cases are correct
          Log::info() << "HybridLinearModelCoeffs<MODEL>::getInfluenceVec() "
                         "unable to get influence region" << std::endl;
          abort();
        }
      }
  }
    return dxVec;
}
//------------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_HYBRIDLINEARMODELCOEFFS_H_
