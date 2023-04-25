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

//------------------------------------------------------------------------------

template <typename MODEL>
class HybridLinearModelCoeffsParameters : public Parameters {
    OOPS_CONCRETE_PARAMETERS(HybridLinearModelCoeffsParameters, Parameters)
    typedef HtlmEnsembleParameters<MODEL>               HtlmEnsembleParameters_;
    typedef HtlmCalculatorParameters<MODEL>             HtlmCalculatorParameters_;
    typedef typename Geometry<MODEL>::Parameters_       GeometryParameters_;
    typedef typename Increment<MODEL>::ReadParameters_  ReadParameters_;
 public:
    // Ensemble parameters
    RequiredParameter<HtlmEnsembleParameters_> htlmEnsemble{"htlm ensemble", this};
    // Coefficent calculator parameters
    RequiredParameter<HtlmCalculatorParameters_> htlmCalculator{"htlm calculator", this};
    // Assimilation window length
    RequiredParameter<util::Duration> windowLength{"window length", this};
    // Window begin
    RequiredParameter<util::DateTime> windowBegin{"window begin", this};
    // Variables
    RequiredParameter<Variables> vars{"training variables", this};
};

//------------------------------------------------------------------------------

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

// The geometry for the state is a .yaml parameter
// This is used for the TLM geometry in HtlmEnsemble.h
// This is ultimately first constructed in HybridLinearModel.h from its parameters
    HybridLinearModelCoeffs(const HybridLinearModelCoeffsParameters_ &,
                            const Geometry_ &,
                            const util::Duration &);

 public:
    void updateIncTL(Increment_ &) const;
    void updateIncAD(Increment_ &) const;

 private:
    atlas::Field makeUpdateStencil(const atlas::idx_t &);

 private:
    std::map<util::DateTime, atlas::FieldSet> coeffSaver_;
    const Variables vars_;
    const atlas::idx_t influenceSize_;
    const atlas::idx_t halfInfluenceSize_;
    atlas::Field       influenceStencil_;
};

//------------------------------------------------------------------------------

template<typename MODEL>
HybridLinearModelCoeffs<MODEL>::HybridLinearModelCoeffs
  (const HybridLinearModelCoeffsParameters_ & params,
   const Geometry_ & geomTLM,
   const util::Duration & tstep) :
  vars_(params.vars.value()),
  influenceSize_(params.htlmCalculator.value().influenceRegionSize),
  halfInfluenceSize_(influenceSize_ / 2) {
    const util::DateTime windowBegin = params.windowBegin.value();
    const util::Duration windowLength = params.windowLength.value();
    util::DateTime time(windowBegin);
    Increment_ vertExtInc(geomTLM , vars_ , time);
    atlas::FieldSet vertExtFset = vertExtInc.fieldSet();
    atlas::idx_t vertExt_ = vertExtFset[0].shape(1);
    influenceStencil_ = makeUpdateStencil(vertExt_);
    HtlmEnsemble_ ens(params.htlmEnsemble.value(), geomTLM);
    HtlmCalculator_ calculator(params.htlmCalculator.value(), params.vars,
    params.htlmEnsemble.value().ensembleSize.value(), geomTLM, params.windowBegin.value());
    const size_t ensembleSize = params.htlmEnsemble.value().ensembleSize.value();
    // Step ensemble and calculate coefficients at time t
    while (time < (windowBegin + windowLength)) {
      time += tstep;
      // Step ensemble (HtlmEnsemble.h)
      ens.step(tstep);
      // Create an empty atlas::FieldSet (coeffFieldSet) to store coefficients at this time
      // Dimensions of the FieldSet are defined in the HtlmCalculator object
      // (To which the field is passed by reference)
      atlas::FieldSet coeffFieldSet;
      // Emplace coeffFieldSet in map for storage
      coeffSaver_.emplace(time, coeffFieldSet);
      // Ensemble info and FieldSet for coefficient storage passed to HtlmCalculator
      calculator.calcCoeffs(ens.getLinearEns(), ens.getLinearErrDe(), coeffFieldSet);
      // Update increments with coefficients before next step
      for (size_t ensInd = 0; ensInd < ensembleSize; ++ensInd) {
          updateIncTL(ens.getLinearEns()[ensInd]);
      }
    }
}

//------------------------------------------------------------------------------

// Updates grid point dx to dx' with htlm coefficents
template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::updateIncTL(Increment_ & dx) const {
    Log::trace() << "HybridLinearModelCoeffs<MODEL::updateIncTL() starting" << std::endl;
    auto stencilView = atlas::array::make_view<int, 2>(influenceStencil_);
    atlas::FieldSet & dxFset = dx.fieldSet();
    auto dxView = atlas::array::make_view<double, 2>(dxFset[vars_[0]]);
        for (atlas::idx_t i = 0; i < dxView.shape(0); ++i) {
            std::vector<double> updateVals(dxView.shape(1)*vars_.size(), 0.0);
            // Calculate update values
            for (size_t varInd = 0; varInd < vars_.size(); ++varInd) {
                auto coeffView = atlas::array::make_view<double, 3>
                                          (coeffSaver_.at(dx.validTime())[vars_[varInd]]);
                for (atlas::idx_t k = 0; k < dxView.shape(1); k++) {
                    for (size_t varInd2 = 0; varInd2 < vars_.size(); varInd2++) {
                        auto dxView = atlas::array::make_view<double, 2>(dxFset[vars_[varInd2]]);
                        for (atlas::idx_t coeffInd = 0; coeffInd < influenceSize_; coeffInd++) {
                            updateVals[k+varInd*dxView.shape(1)] +=
                                coeffView(i, k, varInd2*influenceSize_ + coeffInd)*
                                                      dxView(i, stencilView(k, coeffInd));
                        }
                    }
                }
            }
            // Update column
            for (size_t varInd = 0; varInd < vars_.size(); ++varInd) {
                auto dxView = atlas::array::make_view<double, 2>(dxFset[vars_[varInd]]);
                for (atlas::idx_t k = 0; k < dxView.shape(1); k++) {
                    dxView(i, k)+=updateVals[k+varInd*dxView.shape(1)];
                }
            }
        }
    dx.synchronizeFields();
    Log::trace() << "HybridLinearModelCoeffs<MODEL::updateIncTL() done" << std::endl;
}

//------------------------------------------------------------------------------

// Adjoint of updateIncTL
template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::updateIncAD(Increment_ & dx) const {
    Log::trace() << "HybridLinearModelCoeffs<MODEL::updateIncAD() starting" << std::endl;
    auto stencilView = atlas::array::make_view<int, 2>(influenceStencil_);
    atlas::FieldSet & dxFset = dx.fieldSet();
    auto dxView = atlas::array::make_view<double, 2>(dxFset[vars_[0]]);
        // Adjoint statement for updateVals
        for (atlas::idx_t i = 0; i < dxView.shape(0); ++i) {
            std::vector<double> updateVals(dxView.shape(1)*vars_.size(), 0.0);
            for (size_t varInd = 0; varInd < vars_.size(); ++varInd) {
                auto dxView = atlas::array::make_view<double, 2>(dxFset[vars_[varInd]]);
                for (atlas::idx_t k = 0; k < dxView.shape(1); k++) {
                    updateVals[k+varInd*dxView.shape(1)]+=dxView(i, k);
                }
            }
            // Adjoint statement for dxView
            for (size_t varInd = 0; varInd < vars_.size(); ++varInd) {
                auto coeffView = atlas::array::make_view<double, 3>
                                          (coeffSaver_.at(dx.validTime())[vars_[varInd]]);
                for (atlas::idx_t k = 0; k < dxView.shape(1); k++) {
                    for (size_t varInd2 = 0; varInd2 < vars_.size(); varInd2++) {
                        auto dxView = atlas::array::make_view<double, 2>(dxFset[vars_[varInd2]]);
                        for (atlas::idx_t coeffInd = 0; coeffInd < influenceSize_; coeffInd++) {
                            dxView(i, stencilView(k, coeffInd)) +=
                                    coeffView(i, k, varInd2*influenceSize_ + coeffInd)*
                                                        updateVals[k+varInd*dxView.shape(1)];
                        }
                    }
                }
            }
        }
    dx.synchronizeFields();
    Log::trace() << "HybridLinearModelCoeffs<MODEL::updateIncAD() done" << std::endl;
}

//------------------------------------------------------------------------------

template<typename MODEL>
atlas::Field HybridLinearModelCoeffs<MODEL>::makeUpdateStencil(const atlas::idx_t & vertExt_) {
    atlas::Field
         influenceStencil_("influence_stencil", atlas::array::make_datatype<int>(),
            atlas::array::make_shape(vertExt_, influenceSize_));
    auto stencilView = atlas::array::make_view<atlas::idx_t, 2>(influenceStencil_);
    for (atlas::idx_t k = 0; k < vertExt_; ++k) {
        for (atlas::idx_t infInd = 0; infInd < influenceSize_; ++infInd) {
            if (k-halfInfluenceSize_ > 0 &&
                    k + halfInfluenceSize_ < vertExt_ ) {
            // middle case
            stencilView(k, infInd) = k-halfInfluenceSize_ + infInd;
            } else if (k-halfInfluenceSize_ <= 0) {
            // start of increment edge case
            stencilView(k, infInd) = infInd;
            } else if (k+halfInfluenceSize_ >= vertExt_) {
            // end of increment edge case
            stencilView(k, infInd) = (vertExt_-influenceSize_) + infInd;
            } else {
            // Checking cases are correct
            Log::info() << "HtlmLinearModelCoeffs<MODEL>::makeUpdateStencil() "
                            "unable to get influence region" << std::endl;
            abort();
            }  // end else
        }   // end infInd
    }
return influenceStencil_;
}

}  // namespace oops

#endif  // OOPS_GENERIC_HYBRIDLINEARMODELCOEFFS_H_
