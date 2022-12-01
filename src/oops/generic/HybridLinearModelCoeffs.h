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
    // Precalculated coefficient files (if using)
    Parameter<std::vector<ReadParameters_>> coefficients{"coefficients", {}, this};
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
    std::vector<double> getInfluenceVec(const atlas::FieldSet &,
                                        const atlas::idx_t & i,
                                        const atlas::idx_t & k) const;

 private:
    std::map<util::DateTime, atlas::FieldSet> coeffSaver_;
    const Variables vars_;
    const atlas::idx_t influenceSize_;
    const atlas::idx_t halfInfluenceSize_;
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
    if (!params.coefficients.value().empty()) {
        size_t fileIndex = 0;
        while (time < (windowBegin + windowLength)) {
            time += tstep;
            atlas::FieldSet rank3FieldSet;  // Equivalent to coeffFieldSet in else clause below
            for (atlas::idx_t x = 0; x < influenceSize_; x++) {
                Increment_ incrementFromFile(geomTLM, vars_, time);
                incrementFromFile.read(params.coefficients.value()[fileIndex]);
                ++fileIndex;
                atlas::FieldSet & rank2FieldSet = incrementFromFile.fieldSet();
                for (atlas::FieldSet::iterator iterFieldSet = rank2FieldSet.begin();
                     iterFieldSet != rank2FieldSet.end();
                     ++iterFieldSet) {
                    atlas::Field rank2Field = *iterFieldSet;
                    auto rank2FieldView = atlas::array::make_view<double, 2>(rank2Field);
                    atlas::Field rank3Field(rank2Field.name(),
                                            atlas::array::make_datatype<double>(),
                                            atlas::array::make_shape(rank2Field.shape(0),
                                                                     rank2Field.shape(1),
                                                                     influenceSize_));
                    auto rank3FieldView = atlas::array::make_view<double, 3>(rank3Field);
                    for (atlas::idx_t i = 0; i < rank2Field.shape(0); i++) {
                        for (atlas::idx_t k = 0; k < rank2Field.shape(1); k++) {
                            rank3FieldView(i, k, x) = rank2FieldView(i, k);
                        }
                    }
                    rank3FieldSet.add(rank3Field);
                }
            }
            coeffSaver_.emplace(time, rank3FieldSet);
        }
    } else {
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
}

//------------------------------------------------------------------------------

// Updates grid point dx to dx' via dx' = dx + dotproduct(coefficent vector, dx influence region)
template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::updateIncTL(Increment_ & dx) const {
    Log::trace() << "HybridLinearModelCoeffs<MODEL::updateIncTL() starting" << std::endl;
    atlas::FieldSet & dxFset = dx.fieldSet();
    for (size_t varInd = 0; varInd < vars_.size(); ++varInd) {
        auto dxView = atlas::array::make_view<double, 2>(dxFset[vars_[varInd]]);
        for (atlas::idx_t i = 0; i < dxView.shape(0); ++i) {
            std::vector<double> updateVal(dxView.shape(1), 0.0);
            for (atlas::idx_t k = 0; k < dxView.shape(1); ++k) {
                std::vector<double> dxVec = getInfluenceVec(dxFset, i, k);
                auto coeffView = atlas::array::make_view<double, 3>
                        (coeffSaver_.at(dx.validTime())[vars_[varInd]]);
                for (atlas::idx_t infInd = 0; infInd < atlas::idx_t(dxVec.size()); ++infInd) {
                    updateVal[k] += coeffView(i, k, infInd) * dxVec[infInd];
                }
            }
            for (atlas::idx_t k = 0; k < dxView.shape(1); ++k) {
                dxView(i, k) += updateVal[k];
            }
        }
    }
    dx.synchronizeFields();
    Log::trace() << "HybridLinearModelCoeffs<MODEL::updateIncTL() done" << std::endl;
}

//------------------------------------------------------------------------------

// Adjoint of above
template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::updateIncAD(Increment_ & dx) const {
    Log::trace() << "HybridLinearModelCoeffs<MODEL::updateIncAD() starting" << std::endl;
    atlas::FieldSet & dxFset = dx.fieldSet();
    for (size_t varInd = 0; varInd < vars_.size(); ++varInd) {
        auto dxView = atlas::array::make_view<double, 2>(dxFset[vars_[varInd]]);
        for (atlas::idx_t i = 0; i < dxView.shape(0); ++i) {
            std::vector<double> updateVal(dxView.shape(1), 0.0);
            for (atlas::idx_t k = 0; k < dxView.shape(1); ++k) {
                auto coeffView = atlas::array::make_view<double, 3>
                        (coeffSaver_.at(dx.validTime())[vars_[varInd]]);
                if (k - halfInfluenceSize_ > 0 && k + halfInfluenceSize_ < dxView.shape(1)) {
                    for (int infInd = 0; infInd < influenceSize_; infInd++) {
                        updateVal[k + infInd - halfInfluenceSize_] +=
                                coeffView(i, k, infInd) * dxView(i, k);
                    }
                } else if (k - halfInfluenceSize_ <= 0) {
                    for (int infInd = 0; infInd < influenceSize_; infInd++) {
                        updateVal[infInd] += coeffView(i, k, infInd) * dxView(i, k);
                    }
                } else if (k + halfInfluenceSize_ >= dxView.shape(1)) {
                    for (int infInd= 0; infInd < influenceSize_; infInd++) {
                        updateVal[(dxView.shape(1) - influenceSize_ + infInd)] +=
                                coeffView(i, k, infInd) * dxView(i, k);
                    }
                } else {
                    // Checking cases are correct
                    Log::info() << "HybridLinearModelCoeffs<MODEL>::updateIncAD()"
                                   "unable to get influence region AD" << std::endl;
                    abort();
                }
            }
            for (atlas::idx_t k = 0; k < dxView.shape(1); ++k) {
                dxView(i, k) += updateVal[k];
            }
        }
    }
    dx.synchronizeFields();
    Log::trace() << "HybridLinearModelCoeffs<MODEL::updateIncAD() done" << std::endl;
}

//------------------------------------------------------------------------------

// Create influence region vector for grid point dx (gets points above and below for all variables)
template<typename MODEL>
std::vector<double> HybridLinearModelCoeffs<MODEL>::getInfluenceVec(const atlas::FieldSet & dxFset,
                                                                    const atlas::idx_t & i,
                                                                    const atlas::idx_t & k) const {
    Log::trace() << "HybridLinearModelCoeffs<MODEL>::getInfluenceVec() starting" << std::endl;
    std::vector<double> dxVec(influenceSize_ * vars_.size());
    for (size_t varInd = 0; varInd < vars_.size(); ++varInd) {
        const auto dxFview = atlas::array::make_view<double, 2>(dxFset[vars_[varInd]]);
        for (atlas::idx_t infInd = 0; infInd < influenceSize_; ++infInd) {
            if (k - halfInfluenceSize_ > 0 && k + halfInfluenceSize_ < dxFview.shape(1)) {
                // General (middle) case
                dxVec[vars_.size() * varInd + infInd] =
                        dxFview(i, k - halfInfluenceSize_ + infInd);
            } else if (k - halfInfluenceSize_ <= 0) {
                // Start of increment edge case
                dxVec[vars_.size() * varInd + infInd] = dxFview(i, infInd);
            } else if (k + halfInfluenceSize_ >= dxFview.shape(1)) {
                // End of increment edge case
                dxVec[vars_.size() * varInd + infInd] =
                        dxFview(i, (dxFview.shape(1) - influenceSize_) + infInd);
            } else {
                // Checking cases are correct
                Log::info() << "HybridLinearModelCoeffs<MODEL>::getInfluenceVec()"
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
