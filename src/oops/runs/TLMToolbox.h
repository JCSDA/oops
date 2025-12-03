/*
 * (C) Copyright 2023-2025 UCAR.
 * (C) Crown copyright 2023-2025 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_TLMTOOLBOX_H_
#define OOPS_RUNS_TLMTOOLBOX_H_

#include <optional>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "oops/base/LinearModel.h"
#include "oops/base/Model.h"
#include "oops/base/StructuredGridWriter.h"
#include "oops/base/TrajectorySaver.h"
#include "oops/generic/instantiateLinearModelFactory.h"
#include "oops/generic/instantiateModelFactory.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/dot_product.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/PrintAdjTest.h"
#include "oops/util/printRunStats.h"
#include "oops/util/TimeWindow.h"

namespace oops {

/// \brief Run tangent linear model forecast, optionally compute linearization error, test adjoint.
///
/// \details
///
/// Configuration for a tangent linear model forecast:
/// - "geometry:" OR ("linear geometry:" AND "nonlinear geometry:")   - Same or separate geometries.
/// - "linear model:"
/// - "model aux increment:"                                          - Set to null if not required.
/// - "model:"                                                            - For trajectory forecast.
/// - "model aux control:"                                            - Set to null if not required.
/// - "x1:"                                 - Nonlinear model trajectory forecast initial condition.
/// - "x2:" OR "dx:"                    - Initial increment computed as (x2 - x1) or read from file.
/// - "linear variables:"     - Optional; variables for linear model; taken from x1 if not provided.
/// - "forecast length:"
/// - "time resolution:"          - Optional frequency to save dx; must factor into forecast length.
/// - "output:
///      dx:
///        high res:                      - Optional; default false; save dx at nonlinear geometry*.
///        model grid:                             - Optional; save dx using model interface writer.
///        structured grid:                          - Optional; save dx using StructuredGridWriter.
///      error:               "           - Must include, but set to null since not computing error.
///
/// *Cannot be true if saving dx on a structured grid; would invoke an unnecessary interpolation.
///
/// Configuration for a linearization error computation is the same as above, with these additions:
/// - Must use "x2:", not "dx:", so that there is a nonlinear difference to compute error against
/// - "compute error:" - Set to true.
/// - "output:
///      error:
///        high res:           - Optional; default true; compute/save error at nonlinear geometry*.
///        model grid:
///        structured grid:   "
///
/// *This option is more important in the case of the error, as it defines how the error is computed
/// as well as what grid it is saved on. Details below. Unlike in the case of the forecast, it can
/// be true even when saving error on a structured grid.
///
/// Linearization error is defined as:
///   L * dx - ( N(x + dx) - N(x) )
/// where:
///   L is a linear operator representing a linear model forecast of length T,
///   dx is a small perturbation (e.g., analysis increment) valid at the initial time,
///   N is an operator representing a nonlinear model forecast of length T,
///   x is a model state valid at the initial time.
///
/// L is the linearization of N about the trajectory N(x).
///
/// The equations above assume that both L and N are on the same grid. If this is not true,
/// regridding is required. This produces two options for computing error; on the grid of L, or on
/// the grid of N. Either option is available in this application. The former is:
///   L * P1 * (x2 - x1) - P1 * ( N(x2) - N(x1) )
/// where:
///   P1 is a linear operator representing a transformation from the N grid to the L grid,
///   x2 is a perturbed model state valid at the initial time (i.e., x + dx),
///   x1 is an unperturbed model state valid at the initial time (i.e., x).
/// The latter is:
///   P2 * L * P1 * (x2 - x1) - ( N(x2) - N(x1) )
/// where:
///   P2 is a linear operator representing a transformation from the L grid to the N grid.
/// The second option is the default in this application.
///
/// Configuration for an adjoint test is the same as for a tangent linear model forecast, except:
/// - Exclude "x2:"/"dx:"; this indicates that an adjoint test is desired and the increments will be
///   initialised randomly.
/// - Do not add "time resolution:"; the adjoint test is done for "forecast length:".
/// - Set both "output.dx:" and "output.error:" to null.
///
/// Field normalization is employed to ensure that the adjoint test is fair/avoid false positives:
///   dx <- random
///   Mdx <- M * dx
///   S <- Diagonal matrix of per-variable inner products of Mdx
///   SMdx <- S * Mdx
///   dy <- random
///   Sdy <- S * dy
///   ASdy <- A * Sdy
///   Compare inner products <SMdx, dy> == <dx, ASdy>
/// Note variable names in the code differ slightly from those used here (for more optimal code).
/// There is no pass/fail criterion - inner product comparisons are output to the logs.
template <typename MODEL> class TLMToolbox : public Application {
    typedef Geometry<MODEL>                     Geometry_;
    typedef typename std::shared_ptr<Geometry_> PGeometry_;
    typedef Increment<MODEL>                    Increment_;
    typedef LinearModel<MODEL>                  LinearModel_;
    typedef Model<MODEL>                        Model_;
    typedef ModelAuxControl<MODEL>              ModelAuxControl_;
    typedef ModelAuxIncrement<MODEL>            ModelAuxIncrement_;
    typedef PostProcessorTLAD<MODEL>            PostProcessorTLAD_;
    typedef State<MODEL>                        State_;
    typedef StructuredGridWriter<MODEL>         StructuredGridWriter_;
    typedef TrajectorySaver<MODEL>              TrajectorySaver_;

 public:
    explicit TLMToolbox(const eckit::mpi::Comm& comm = oops::mpi::world()) : Application(comm) {
        instantiateLinearModelFactory<MODEL>();
        instantiateModelFactory<MODEL>();
    }

    virtual ~TLMToolbox() = default;

// -------------------------------------------------------------------------------------------------

    int execute(const eckit::Configuration& config) const override {
        Log::trace() << "TLMToolbox::execute() start" << std::endl;
        util::printRunStats("TLMToolbox start");

        // Validate configuration and determine if linearization error is to be computed
        validateConfiguration(config);
        const bool computeError = config.getBool("compute error", false);
        const eckit::LocalConfiguration linearModelConfig(config, "linear model");
        const eckit::LocalConfiguration modelConfig(config, "model");
        const eckit::LocalConfiguration x1Config(config, "x1");

        // Set up geometry(ies)
        const auto[linearGeometry, nonlinearGeometry] = createGeometries(config);

        // Set up linear model (shared pointer must be used for TrajectorySaver constructor)
        const std::shared_ptr<LinearModel_> linearModel = std::make_shared<LinearModel_>(
            *linearGeometry, linearModelConfig);
        ModelAuxIncrement_ mAuxInc(
            *linearGeometry, config.getSubConfiguration("model aux increment"));

        // Set up nonlinear model
        const Model_ model(*nonlinearGeometry, modelConfig);
        const ModelAuxControl_ mAuxCtl(
            *nonlinearGeometry, config.getSubConfiguration("model aux control"));

        // Set up trajectory saver for use with x1 (and empty post-processor for x2 if required)
        const PostProcessorTLAD_ ppT;
        auto[trajSaver, emptyPp] = createPostProcessors(
            linearModelConfig, *linearGeometry, mAuxCtl, linearModel, ppT);

        // Set up trajectory initial condition x1
        State_ x1(*nonlinearGeometry, x1Config);

        // Set up perturbed initial condition x2 if required
        std::optional<State_> x2;
        if (config.has("x2")) x2 = State_(*nonlinearGeometry, config.getSubConfiguration("x2"));

        // Set up increment dx, either as (x2 - x1), from file, or random (for adjoint test only)
        const Variables linearVariables(config.getStringVector("linear variables",
                                                               x1.variables().variables()));
        Increment_ dxHighRes(*nonlinearGeometry, linearVariables, x1.validTime());
        Increment_ dx = createDx(config, *linearGeometry, linearVariables, x1, *x2, dxHighRes);

        Log::test() << "dx at " << dx.validTime() << ":" << dx << std::endl;

        // If running adjoint test, save initial dx
        std::optional<Increment_> dxInitial;
        if (!config.has("x2") && !config.has("dx")) dxInitial = Increment_(dx);

        // Set up time, forecast length and time resolution at which to (compute linearization
        // error and) save fields
        util::DateTime time(dx.validTime());
        const util::Duration forecastLength(config.getString("forecast length"));
        const util::DateTime endTime = time + forecastLength;
        const util::Duration timeResolution(
            config.getString("time resolution", forecastLength.toString()));

        // Set up I/O
        const bool saveDxAtHighRes
            = config.getBool("output.dx.high res", false);
        const bool saveErrorAtHighRes
            = config.getBool("output.error.high res", true);
        const Writer dxWriter(config.getSubConfiguration("output.dx"), *linearGeometry);
        const Writer errorWriter(
            config.getSubConfiguration("output.error"),
            saveErrorAtHighRes ? *nonlinearGeometry : *linearGeometry);

        // Loop over steps of length timeResolution until window is complete
        while (time < endTime) {
            // Forecast x1 from time to time + timeResolution using model
            model.forecast(x1, mAuxCtl, timeResolution, trajSaver);

            // Forecast dx from time to time + timeResolution using linear model
            linearModel->forecastTL(dx, mAuxInc, timeResolution);

            time += timeResolution;
            Log::test() << "dx at " << time << ":" << dx << std::endl;

            // Save
            if (saveDxAtHighRes) {
                dxHighRes = Increment_(*nonlinearGeometry, dx);
                dxWriter.write(dxHighRes);
            } else {
                dxWriter.write(dx);
            }

            if (computeError) {
                // Forecast x2 from time to time + timeResolution using model
                model.forecast(*x2, mAuxCtl, timeResolution, emptyPp);

                // Compute difference between two states at time + timeResolution
                dxHighRes.updateTime(timeResolution);
                dxHighRes.diff(*x2, x1);

                // Compute and save error
                if (saveErrorAtHighRes) {
                    Increment_ error(*nonlinearGeometry, dx);
                    error -= dxHighRes;
                    Log::test() << "error at " << time << ":" << error << std::endl;
                    errorWriter.write(error);
                } else {
                    Increment_ error(dx);
                    error -= Increment_(*linearGeometry, dxHighRes);
                    Log::test() << "error at " << time << ":" << error << std::endl;
                    errorWriter.write(error);
                }
            }
        }

        // Adjoint test
        if (!config.has("x2") && !config.has("dx")) {
            // Normalize each variable of dx using per-variable inner products
            const auto innerProducts = innerProductByVariable(dx);
            auto dxNormalized(dx);
            normalize(dxNormalized, innerProducts);
            // Initialise random increment dy
            auto dy(dx);
            dy.random();
            // Scale dy using normalization factors from dx
            auto dyNormalized(dy);
            normalize(dyNormalized, innerProducts);
            // Adjoint forecast
            linearModel->forecastAD(dyNormalized, mAuxInc, forecastLength);
            // Compare <SMdx, dy> == <dx, ASdy> (S is diagonal matrix of normalization factors)
            const auto inner1 = dot_product(dxNormalized, dy);
            const auto inner2 = dot_product(*dxInitial, dyNormalized);
            Log::info() << util::PrintAdjTest(inner1, inner2, "SM");
            Log::test() << util::PrintAdjTest(inner1, inner2, "SM");
        }

        util::printRunStats("TLMToolbox end");
        Log::trace() << "TLMToolbox::execute() done" << std::endl;
        return 0;
    }

// -------------------------------------------------------------------------------------------------

 private:
    std::string appname() const override { return "oops::TLMToolbox<" + MODEL::name() + ">"; }

// -------------------------------------------------------------------------------------------------

    void validateConfiguration(const eckit::Configuration& c) const {
        if (!((c.has("geometry") && !c.has("linear geometry") && !c.has("nonlinear geometry")) ||
              (c.has("linear geometry") && c.has("nonlinear geometry") && !c.has("geometry")))) {
            throw eckit::BadParameter("TLMToolbox: define either \"geometry\" or "
                                      "\"linear geometry\" and \"nonlinear geometry\".");
        }

        if (!((c.has("x2") && !c.has("dx")) ||
              (c.has("dx") && !c.has("x2")) ||
              (!c.has("dx") && !c.has("x2")))) {
            throw eckit::BadParameter("TLMToolbox: define either \"x2\", \"dx\" or neither.");
        }

        if (!c.has("forecast length")) {
            throw eckit::BadParameter("TLMToolbox: define \"forecast length\".");
        }

        if (c.has("time resolution")) {
            if (!(c.has("x2") || c.has("dx"))) {
                throw eckit::BadParameter("TLMToolbox: if defining \"time resolution\", "
                                          "define either \"x2\", \"dx\".");
            }
            if (!(util::Duration(c.getString("forecast length"))
                  % util::Duration(c.getString("time resolution")) == 0)) {
                throw eckit::BadParameter("TLMToolbox: \"forecast length\" not divisible by "
                                          "\"time resolution\".");
            }
        }

        if (c.getBool("compute error", false) && !c.has("x2")) {
            throw eckit::BadParameter("TLMToolbox: if \"compute error\" is true, define \"x2\".");
        }

        if (c.getBool("output.dx.high res", false)) {
            if (c.has("output.dx.structured grid")) {
                throw eckit::BadParameter("TLMToolbox: cannot have \"output.dx.structured grid\" "
                                          "if \"output.dx.high res\" is true.");
            }
        }
    }

// -------------------------------------------------------------------------------------------------

    std::pair<PGeometry_, PGeometry_> createGeometries(const eckit::Configuration& config) const {
        // The geometries can either be the same or different. Use shared pointers so that two
        // identical objects don't have to be constructed if they're the same.
        PGeometry_ linearGeometry, nonlinearGeometry;
        if (config.has("geometry")) {
            linearGeometry = std::make_shared<Geometry_>(
                config.getSubConfiguration("geometry"), this->getComm());
            nonlinearGeometry = linearGeometry;
        } else {
            linearGeometry = std::make_shared<Geometry_>(
                config.getSubConfiguration("linear geometry"), this->getComm());
            nonlinearGeometry = std::make_shared<Geometry_>(
                config.getSubConfiguration("nonlinear geometry"), this->getComm());
        }
        return {linearGeometry, nonlinearGeometry};
    }

// -------------------------------------------------------------------------------------------------

    std::pair<PostProcessor<State_>, PostProcessor<State_>> createPostProcessors(
        const eckit::Configuration& linearModelConfig, const Geometry_& geometry,
        const ModelAuxControl_& mAuxCtl, std::shared_ptr<LinearModel_> linearModel,
        PostProcessorTLAD_ ppT) const {
            PostProcessor<State_> trajSaver;
            trajSaver.enrollProcessor(new TrajectorySaver_(
                linearModelConfig, geometry, mAuxCtl, linearModel, ppT));
            PostProcessor<State_> emptyPp;
            return {trajSaver, emptyPp};
    }

// -------------------------------------------------------------------------------------------------

    Increment_ createDx(const eckit::Configuration& config, const Geometry_& linearGeometry,
        const Variables& linearVariables, const State_& x1, const State_& x2,
        Increment_& dxHighRes) const {
            // If user has defined a perturbed nonlinear model forecast, the initial increment is
            // the difference between the perturbed and unperturbed initial conditions
            if (config.has("x2")) {
                dxHighRes.diff(x2, x1);
                return Increment_(linearGeometry, dxHighRes);
            // If the user has defined an initial increment explicitly, use that
            } else if (config.has("dx")) {
                const eckit::LocalConfiguration dxConfig(config, "dx");
                // It can be read in at either the nonlinear or linear model resolution
                if (dxConfig.getBool("high res", false)) {
                    dxHighRes.read(dxConfig);
                    return Increment_(linearGeometry, dxHighRes);
                } else {
                    Increment_ dx(linearGeometry, linearVariables, x1.validTime());
                    dx.read(dxConfig);
                    return dx;
                }
            // If neither a perturbed nonlinear model forecast nor an initial increment was defined,
            // the user is requesting an adjoint test, so the initial increment is randomised
            } else {
                Increment_ dx(linearGeometry, linearVariables, x1.validTime());
                dx.random();
                return dx;
            }
    }

// -------------------------------------------------------------------------------------------------

    class Writer {
     public:
        Writer(const eckit::LocalConfiguration& config, const Geometry_& geometry)
        : writer_(config.has("structured grid") ?
              std::make_optional<StructuredGridWriter_>(
                  config.getSubConfiguration("structured grid"), geometry)
              : std::nullopt),
          writeConfig_(config.getSubConfiguration("model grid")) {}

        void write(const Increment_& inc) const {
            if (writer_) writer_->interpolateAndWrite(inc);
            if (!writeConfig_.empty()) inc.write(writeConfig_);
        }

     private:
        const std::optional<StructuredGridWriter_> writer_;
        const eckit::LocalConfiguration writeConfig_;
    };

// -------------------------------------------------------------------------------------------------

    std::unordered_map<std::string, double> innerProductByVariable(const Increment_& dx) const {
        if (dx.geometry().fields().empty()) {
            ABORT("TLMToolbox: adjoint test requires Atlas interface (for normalization)");
        }

        std::vector<double> innerProductsVector;
        innerProductsVector.reserve(dx.fieldSet().size());
        for (const auto& field : dx.fieldSet()) {
            // TODO(tom-j-h) get dotProductFields to ignore masked points
            if (field.metadata().has("mask")) {
                ABORT("TLMToolbox: " + field.name() + " has mask, currently not supported");
            }
            innerProductsVector.push_back(util::dotProductFields(field, field,
                                                                 dx.geometry().getComm()));
        }
        dx.geometry().getComm().allReduceInPlace(
            innerProductsVector.begin(), innerProductsVector.end(), eckit::mpi::sum());

        std::unordered_map<std::string, double> innerProducts;
        size_t f = 0;
        for (const auto& field : dx.fieldSet()) {  // repeat iteration over FieldSet to retain order
            innerProducts[field.name()] = innerProductsVector[f];
            f++;
        }

        return innerProducts;
    }

// -------------------------------------------------------------------------------------------------

    void normalize(Increment_& dx, const std::unordered_map<std::string, double>& factors) const {
        for (auto& field : dx.fieldSet()) {
            if (factors.at(field.name()) == 0) {
                ABORT("TLMToolbox: normalization factor for " + field.name() + " is zero");
            }
            const auto scalingFactor = (1.0 / factors.at(field.name()));
            util::multiplyField(field, scalingFactor);
        }
        dx.synchronizeFields();
    }

// -------------------------------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_TLMTOOLBOX_H_
