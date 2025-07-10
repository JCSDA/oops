/*
 * (C) Copyright 2025-2025 UCAR
 * (C) Copyright 2025-2025 NOAA/EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Geometry.h"
#include "oops/base/State.h"
#include "oops/base/StructuredGridWriter.h"
#include "oops/generic/GlobalInterpolator.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"

namespace oops {

// -------------------------------------------------------------------------------------------------
/**
 * @brief Application for interpolating states between two different models
 *
 * This class provides functionality to merge states from different models. The target's original
 * state fields are preserved, and the source's fields are interpolated onto the target's grid and
 * added to the target's state. The application allows either model to be the source or target of
 * the interpolation. The interpolation uses OOPS's unstructured grid interpolator.
 */
template <typename MODEL1, typename MODEL2> class InterpolateStateBetweenModels :
                                                  public Application {
  typedef Geometry<MODEL1>      Geometry1_;
  typedef Geometry<MODEL2>      Geometry2_;
  typedef State<MODEL1>         State1_;
  typedef State<MODEL2>         State2_;

 public:
// -------------------------------------------------------------------------------------------------
  explicit InterpolateStateBetweenModels(const eckit::mpi::Comm & comm = oops::mpi::world())
           : Application(comm) {}
// -------------------------------------------------------------------------------------------------
  virtual ~InterpolateStateBetweenModels() = default;
// -------------------------------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
    // Initialize configurations and models
    const eckit::LocalConfiguration conf1(fullConfig, MODEL1::name());
    const eckit::LocalConfiguration conf2(fullConfig, MODEL2::name());

    // Setup the geometries
    const eckit::LocalConfiguration geom1Config(conf1, "geometry");
    const Geometry1_ geom1(geom1Config, this->getComm());
    const eckit::LocalConfiguration geom2Config(conf2, "geometry");
    const Geometry2_ geom2(geom2Config, this->getComm());

    // Read in the states
    const eckit::LocalConfiguration state1Config(conf1, "state");
    State1_ xx1(geom1, state1Config);
    oops::Log::test() << MODEL1::name() << " state: " << xx1 << std::endl;

    const eckit::LocalConfiguration state2Config(conf2, "state");
    State2_ xx2(geom2, state2Config);
    oops::Log::test() << MODEL2::name() << " state: " << xx2 << std::endl;

    // For debugging: write original states to lat-lon grids
    if (conf1.has("latlon output"))
      writeLatlonOutput(conf1, geom1, xx1, "latlon output");
    if (conf2.has("latlon output"))
      writeLatlonOutput(conf2, geom2, xx2, "latlon output");

    // Configure interpolator
    eckit::LocalConfiguration conf;
    conf.set("local interpolator type", "oops unstructured grid interpolator");
    const std::string interpSource = fullConfig.getString("interpolator source");

    // Execute the appropriate interpolation based on source
    if (interpSource == MODEL1::name()) {
      interpolateAndSave<MODEL1, MODEL2>(conf, xx1, xx2, conf2);
    } else if (interpSource == MODEL2::name()) {
      interpolateAndSave<MODEL2, MODEL1>(conf, xx2, xx1, conf1);
    } else {
      throw eckit::BadParameter("InterpolateStateBetweenModels: interpolator source must be " +
                               MODEL1::name() + " or " + MODEL2::name(), Here());
    }

    return 0;
  }

 private:
  // Helper method to write state to lat-lon grid
  template <typename MODEL>
  void writeLatlonOutput(const eckit::LocalConfiguration & conf,
                         const Geometry<MODEL> & geom,
                         const State<MODEL> & state,
                         const std::string & configKey) const {
    const eckit::LocalConfiguration latlonConf(conf, configKey);
    const oops::StructuredGridWriter<MODEL> latlon(latlonConf, geom);
    latlon.interpolateAndWrite(state);
  }

  // Helper method to handle the common interpolation and output workflow
  template <typename SOURCE_MODEL, typename TARGET_MODEL>
  void interpolateAndSave(const eckit::LocalConfiguration & interpConf,
                          const State<SOURCE_MODEL> & sourceState,
                          const State<TARGET_MODEL> & targetState,
                          const eckit::LocalConfiguration & targetConf) const {
    // Create and configure interpolator
    oops::GlobalInterpolator interp(interpConf, sourceState.geometry().generic(),
                                    targetState.geometry().functionSpace(),
                                    targetState.geometry().getComm());

    // Prepare fieldsets for interpolation
    atlas::FieldSet xin = targetState.fieldSet().fieldSet();
    atlas::FieldSet xout;

    // Apply interpolation
    interp.apply(sourceState.fieldSet().fieldSet(), xout);

    // Combine variables and fields
    oops::Variables vars = targetState.variables();
    vars += sourceState.variables();
    for (const auto & f : xin) {
      xout.add(f);
    }

    // Create output state
    // Note: xout fieldset at this point has all the fields from both models, with
    // correct vertical levels. The contents of the output state xx_out depend on
    // the TARGET_MODEL implementation of State::fromFieldSet.
    State<TARGET_MODEL> xx_out(targetState.geometry(), vars, targetState.validTime());
    xx_out.fieldSet().fieldSet() = xout;
    xx_out.synchronizeFields();
    oops::Log::test() << TARGET_MODEL::name() << " state after interpolation: " <<
                         xx_out << std::endl;

    // Write interpolated output to lat-lon grid for debugging
    if (targetConf.has("latlon output after interpolation"))
      writeLatlonOutput(targetConf, targetState.geometry(), xx_out,
                        "latlon output after interpolation");

    // Write final state to file
    const eckit::LocalConfiguration outputConfig(targetConf, "output state");
    xx_out.write(outputConfig);
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::InterpolateStateBetweenModels<" + MODEL1::name() + "," + MODEL2::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
