/*
 * (C) Copyright 2022 MetOffice.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_HTLMCOEFFICIENTS_H_
#define OOPS_RUNS_HTLMCOEFFICIENTS_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/Variables.h"
#include "oops/generic/HtlmCalculator.h"
#include "oops/generic/HtlmEnsemble.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL> class HTLMCoefficientsParameters : public ApplicationParameters {
    OOPS_CONCRETE_PARAMETERS(HTLMCoefficientsParameters, ApplicationParameters)
    typedef HtlmEnsembleParameters<MODEL>               HtlmEnsembleParameters_;
    typedef HtlmCalculatorParameters<MODEL>             HtlmCalculatorParameters_;
    typedef typename Geometry<MODEL>::Parameters_       GeometryParameters_;
    typedef typename Increment<MODEL>::WriteParameters_ IncrementWriteParameters_;
 public:
    // Linear model geometry
    RequiredParameter<GeometryParameters_> incrementGeometry{"increment geometry", this};
    // Linear model time step
    RequiredParameter<util::Duration> tstep{"tstep", this};
    // Variables
    RequiredParameter<Variables> vars{"training variables", this};
    // Window length
    RequiredParameter<util::Duration> windowLength{"window length", this};
    // Window begin
    RequiredParameter<util::DateTime> windowBegin{"window begin", this};
    // Ensemble parameters
    RequiredParameter<HtlmEnsembleParameters_> htlmEnsemble{"htlm ensemble", this};
    // Coefficient calculator parameters
    RequiredParameter<HtlmCalculatorParameters_> htlmCalculator{"htlm calculator", this};
    // Increment output parameters
    RequiredParameter<IncrementWriteParameters_> output{"output", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class HTLMCoefficients : public Application {
    typedef HtlmEnsemble<MODEL>                         HtlmEnsemble_;
    typedef HtlmEnsembleParameters<MODEL>               HtlmEnsembleParameters_;
    typedef HtlmCalculator<MODEL>                       HtlmCalculator_;
    typedef HtlmCalculatorParameters<MODEL>             HtlmCalculatorParameters_;
    typedef Geometry<MODEL>                             Geometry_;
    typedef Increment<MODEL>                            Increment_;
    typedef typename Increment<MODEL>::WriteParameters_ IncrementWriteParameters_;
    typedef HTLMCoefficientsParameters<MODEL>           HTLMCoefficientsParameters_;

// -----------------------------------------------------------------------------

 public:
    explicit HTLMCoefficients(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {
    }

// -----------------------------------------------------------------------------

    virtual ~HTLMCoefficients() {}

// -----------------------------------------------------------------------------

    int execute(const eckit::Configuration & fullConfig, bool validate) const override {
        // Calculate coefficients:
        // Deserialize parameters
        HTLMCoefficientsParameters_ params;
        if (validate) params.validate(fullConfig);
        params.deserialize(fullConfig);
        // Set up linear model geometry and time step
        const Geometry_ geomTLM(params.incrementGeometry.value(), this->getComm());
        const util::Duration tstep(params.tstep.value());
        // Set up influence region
        const atlas::idx_t influenceSize_ = (params.htlmCalculator.value().influenceRegionSize);
        const atlas::idx_t halfInfluenceSize_ = (influenceSize_ / 2);
        // Initialise variables
        const Variables vars_ = (params.vars.value());
        // Set up ensemble
        HtlmEnsemble_ ens(params.htlmEnsemble.value(), geomTLM);
        const size_t ensembleSize = params.htlmEnsemble.value().ensembleSize.value();
        // Set up calculator
        HtlmCalculator_ calculator(params.htlmCalculator.value(), params.vars,
                                   params.htlmEnsemble.value().ensembleSize.value(),
                                   geomTLM, params.windowBegin.value());
        // Set up map of util::DateTimes (keys) to rank 3 atlas::FieldSets of coefficients (values)
        std::map<util::DateTime, atlas::FieldSet> coeffSaver_;
        // Set up window and initialise time
        const util::DateTime windowBegin = params.windowBegin.value();
        const util::Duration windowLength = params.windowLength.value();
        util::DateTime time(windowBegin);
        // Loop over time, stepping ensemble forward and calculating coefficients
        while (time < (windowBegin + windowLength)) {
            time += tstep;
            // Step ensemble
            ens.step(tstep);
            // Create an empty atlas::FieldSet to store coefficients at this time
            // Dimensions of the set are defined in the calculator object,
            // to which coeffFieldSet is passed by reference
            atlas::FieldSet coeffFieldSet;
            // Emplace coeffFieldSet in map for storage
            coeffSaver_.emplace(time, coeffFieldSet);
            // Calculate coefficients using ensemble information and store them in coeffFieldSet
            calculator.calcCoeffs(ens.getLinearEns(), ens.getLinearErrDe(), coeffFieldSet);
            // Update increments with coefficients before next step
            for (size_t ensInd = 0; ensInd < ensembleSize; ++ensInd) {
                updateIncTL(ens.getLinearEns()[ensInd],
                            vars_,
                            influenceSize_,
                            halfInfluenceSize_,
                            coeffSaver_);
            }
        }

        // Write to files:
        IncrementWriteParameters_ outputConfig = params.output.value();
        // Loop over objects in map
        for (std::map<util::DateTime, atlas::FieldSet>::iterator iterMap = coeffSaver_.begin();
             iterMap != coeffSaver_.end();
             ++iterMap) {
            util::DateTime dateTime = iterMap->first;
            atlas::FieldSet rank3FieldSet = iterMap->second;
            // Loop over variables
            for (atlas::idx_t x = 0; x < rank3FieldSet.field(0).shape(2); x++) {
                Increment_ saverIncrement(geomTLM, vars_, dateTime);
                atlas::FieldSet & rank2FieldSet = saverIncrement.fieldSet();
                outputConfig.setMember(x + 1);
                // Loop over atlas::FieldSets
                for (atlas::FieldSet::iterator iterFieldSet = rank3FieldSet.begin();
                     iterFieldSet != rank3FieldSet.end();
                     ++iterFieldSet) {
                    atlas::Field rank3Field = *iterFieldSet;
                    auto rank3FieldView = atlas::array::make_view<double, 3>(rank3Field);
                    atlas::Field rank2Field(rank3Field.name(),
                                            atlas::array::make_datatype<double>(),
                                            atlas::array::make_shape(rank3Field.shape(0),
                                                                     rank3Field.shape(1)));
                    auto rank2FieldView = atlas::array::make_view<double, 2>(rank2Field);
                    // Loop over horizontal and vertical coordinates
                    for (atlas::idx_t i = 0; i < rank3Field.shape(0); i++) {
                        for (atlas::idx_t k = 0; k < rank3Field.shape(1); k++) {
                            rank2FieldView(i, k) = rank3FieldView(i, k, x);
                        }
                    }
                    rank2FieldSet.add(rank2Field);
                }
                saverIncrement.synchronizeFields();
                saverIncrement.write(outputConfig);
            }
        }
        return 0;
    }

// -----------------------------------------------------------------------------

    void outputSchema(const std::string & outputPath) const override {
        HTLMCoefficientsParameters_ params;
        params.outputSchema(outputPath);
    }

// -----------------------------------------------------------------------------

    void validateConfig(const eckit::Configuration & fullConfig) const override {
        HTLMCoefficientsParameters_ params;
        params.validate(fullConfig);
    }

// -----------------------------------------------------------------------------

 private:
    std::string appname() const override {
        return "oops::HTLMCoefficients<" + MODEL::name() + ">";
    }

// -----------------------------------------------------------------------------

    void updateIncTL(Increment_ & dx,                     const Variables & vars_,
                     const atlas::idx_t & influenceSize_, const atlas::idx_t & halfInfluenceSize_,
                     std::map<util::DateTime, atlas::FieldSet> & coeffSaver_) const {
        Log::trace() << "updateIncTL() starting" << std::endl;
        atlas::FieldSet & dxFset = dx.fieldSet();
        for (size_t varInd = 0; varInd < vars_.size(); ++varInd) {
            auto dxView = atlas::array::make_view<double, 2>(dxFset[vars_[varInd]]);
            for (atlas::idx_t i = 0; i < dxView.shape(0); ++i) {
                std::vector<double> updateVal(dxView.shape(1), 0.0);
                for (atlas::idx_t k = 0; k < dxView.shape(1); ++k) {
                    std::vector<double> dxVec = getInfluenceVec(dxFset,
                                                                i,
                                                                k,
                                                                vars_,
                                                                influenceSize_,
                                                                halfInfluenceSize_);
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
        Log::trace() << "updateIncTL() done" << std::endl;
    }

// -----------------------------------------------------------------------------

    std::vector<double> getInfluenceVec(const atlas::FieldSet & dxFset,
                                        const atlas::idx_t & i,
                                        const atlas::idx_t & k,
                                        const Variables & vars_,
                                        const atlas::idx_t & influenceSize_,
                                        const atlas::idx_t & halfInfluenceSize_) const {
        Log::trace() << "getInfluenceVec() starting" << std::endl;
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
                    Log::info() << "getInfluenceVec() unable to get influence region" << std::endl;
                    abort();
                }
            }
        }
        return dxVec;
    }
// -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_HTLMCOEFFICIENTS_H_
