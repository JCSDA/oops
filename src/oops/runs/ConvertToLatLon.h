/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_CONVERTTOLATLON_H_
#define OOPS_RUNS_CONVERTTOLATLON_H_

#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/LatLonGridWriter.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/State.h"
#include "oops/base/StateEnsemble.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL> class StateToLatLonParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(StateToLatLonParameters, ApplicationParameters)
  typedef State<MODEL>                   State_;
  typedef Geometry<MODEL>                Geometry_;
  typedef LatLonGridWriterParameters LatLonGridWriterParameters_;

 public:
  typedef typename State_::Parameters_      StateParameters_;
  typedef typename Geometry_::Parameters_   GeometryParameters_;

  RequiredParameter<GeometryParameters_>         stateGeometry{"state geometry", this};
  RequiredParameter<StateParameters_>            state{"state", this};
  RequiredParameter<LatLonGridWriterParameters_> latLonInterp{"latlon interpolation", this};
};

template <typename MODEL> class IncToLatLonParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(IncToLatLonParameters, ApplicationParameters)

  typedef Increment<MODEL>               Increment_;
  typedef Geometry<MODEL>                Geometry_;
  typedef LatLonGridWriterParameters LatLonGridWriterParameters_;

 public:
  typedef typename Increment_::ReadParameters_      IncrementParameters_;
  typedef typename Geometry_::Parameters_           GeometryParameters_;

  RequiredParameter<GeometryParameters_>          incGeometry{"increment geometry", this};
  RequiredParameter<Variables>                    vars{"variables", this};
  RequiredParameter<util::DateTime>               date{"date", this};
  RequiredParameter<IncrementParameters_>         increment{"increment", this};
  RequiredParameter<LatLonGridWriterParameters_>  latLonInterp{"latlon interpolation", this};
};

template <typename MODEL> class StateEnsToLatLonParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(StateEnsToLatLonParameters, ApplicationParameters)

  typedef StateEnsembleParameters<MODEL> StateEnsembleParameters_;
  typedef LatLonGridWriterParameters LatLonGridWriterParameters_;
  typedef Geometry<MODEL>                Geometry_;
 public:
  typedef typename Geometry_::Parameters_ GeometryParameters_;

  RequiredParameter<StateEnsembleParameters_>    stateEnsemble{"states", this};
  RequiredParameter<LatLonGridWriterParameters_> latLonInterp{"latlon interpolation", this};
  RequiredParameter<GeometryParameters_>         stateGeometry{"state geometry", this};
};

template <typename MODEL> class IncEnsToLatLonParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(IncEnsToLatLonParameters, ApplicationParameters)

  typedef IncrementEnsembleParameters<MODEL> IncrementEnsembleParameters_;
  typedef LatLonGridWriterParameters LatLonGridWriterParameters_;
  typedef Geometry<MODEL>                Geometry_;
 public:
  typedef typename Geometry_::Parameters_ GeometryParameters_;

  RequiredParameter<IncrementEnsembleParameters_> incrementEnsemble{"increments", this};
  RequiredParameter<LatLonGridWriterParameters_>  latLonInterp{"latlon interpolation", this};
  RequiredParameter<GeometryParameters_>          incrementGeometry{"increment geometry", this};
  RequiredParameter<Variables>                    incrementVariables{"increment variables", this};
};

template <typename MODEL> class ConvertToLatLonParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(ConvertToLatLonParameters, ApplicationParameters)

  typedef StateEnsToLatLonParameters<MODEL>       StateEnsToLatLonParameters_;
  typedef IncEnsToLatLonParameters<MODEL>         IncEnsToLatLonParameters_;
  typedef StateToLatLonParameters<MODEL>          StateToLatLonParameters_;
  typedef IncToLatLonParameters<MODEL>            IncToLatLonParameters_;

 public:
  OptionalParameter<StateEnsToLatLonParameters_>
                   stateEnsToLatLon{"state ensemble to latlon", this};
  OptionalParameter<IncEnsToLatLonParameters_>
                   incEnsToLatLon{"increment ensemble to latlon", this};
  OptionalParameter<std::vector<StateToLatLonParameters_>>
                   stateToLatLon{"states to latlon", this};
  OptionalParameter<std::vector<IncToLatLonParameters_>>
                   incToLatLon{"increments to latlon", this};
};

template <typename MODEL> class ConvertToLatLon : public Application {
  typedef Geometry<MODEL>                           Geometry_;
  typedef State<MODEL>                              State_;
  typedef StateEnsemble<MODEL>                      StateEnsemble_;
  typedef Increment<MODEL>                          Increment_;
  typedef IncrementEnsemble<MODEL>                  IncrementEnsemble_;
  typedef std::shared_ptr<IncrementEnsemble<MODEL>> EnsemblePtr_;
  typedef LatLonGridWriter<MODEL>                   LatLonGridWriter_;

  typedef ConvertToLatLonParameters<MODEL>          ConvertToLatLonParameters_;

 public:
// -----------------------------------------------------------------------------
  explicit ConvertToLatLon(const eckit::mpi::Comm & comm = oops::mpi::world()) :
                                                              Application(comm) {}
// -----------------------------------------------------------------------------
  virtual ~ConvertToLatLon() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
//  Deserialize parameters
    ConvertToLatLonParameters_ params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

// -----------------------------------------------------------------------------

//  Interpolate state ensemble if provided
    if (params.stateEnsToLatLon.value() != boost::none) {
      Log::info() << "Interpolating State Ensemble" << std::endl;
      Geometry_ resol_(params.stateEnsToLatLon.value()->stateGeometry, this->getComm());
      StateEnsemble_ statesToInterp_(resol_, params.stateEnsToLatLon.value()->
                                                                 stateEnsemble.value());
      const LatLonGridWriter_ latLonWriter_(params.stateEnsToLatLon.value()->latLonInterp.value(),
                                                                                      resol_);
      size_t numstates = statesToInterp_.size();
      for (size_t jm=0; jm < numstates; jm++) {
        latLonWriter_.interpolateAndWrite(statesToInterp_[jm]);
        Log::test() << latLonWriter_ << std::endl;
      }
    }

// -----------------------------------------------------------------------------

//  Interpolate individual states if provided
    if (params.stateToLatLon.value() != boost::none) {
      Log::info() << "Interpolating Individual State(s) " << std::endl;
      size_t numstates = params.stateToLatLon.value()->size();
      for (size_t jm=0; jm < numstates; jm++) {
        Geometry_ resol_(params.stateToLatLon.value()->at(jm).stateGeometry, this->getComm());
        State_ stateToInterp_(resol_, params.stateToLatLon.value()->at(jm).state.value());
        const LatLonGridWriter_ latLonWriter_(params.stateToLatLon.value()->
                                                          at(jm).latLonInterp.value(), resol_);
        latLonWriter_.interpolateAndWrite(stateToInterp_);
        Log::test() << latLonWriter_ << std::endl;
      }
    }

// -----------------------------------------------------------------------------

//  Interpolate individual increments if provided
    if (params.incToLatLon.value() != boost::none) {
      Log::info() << "Interpolating Individual Increment(s)" << std::endl;
      size_t numstates = params.incToLatLon.value()->size();
      for (size_t jm=0; jm < numstates; jm++) {
        Geometry_ resol_(params.incToLatLon.value()->at(jm).incGeometry, this->getComm());
        Increment_ incToInterp_(resol_, params.incToLatLon.value()->at(jm).vars.value(),
                                        params.incToLatLon.value()->at(jm).date.value());
        incToInterp_.read(params.incToLatLon.value()->at(jm).increment.value());
        const LatLonGridWriter_ latLonWriter_(params.stateToLatLon.value()->
                                               at(jm).latLonInterp.value(), resol_);
        // This supports output on model levels only; to output on pressure levels would need to
        // read in a reference background from which to read the vertical pressure coordinate
        latLonWriter_.interpolateAndWrite(incToInterp_);
        Log::test() << latLonWriter_ << std::endl;
      }
    }
// -----------------------------------------------------------------------------

    if (params.incEnsToLatLon.value() != boost::none) {
      Log::info() << "Interpolating Increment Ensemble" << std::endl;
      Geometry_ resol_(params.incEnsToLatLon.value()->incrementGeometry, this->getComm());
      IncrementEnsemble_ incrementsToInterp_(resol_,
                            params.incEnsToLatLon.value()->incrementVariables.value(),
                            params.incEnsToLatLon.value()->incrementEnsemble.value());
      const LatLonGridWriter_ latLonWriter_(params.incEnsToLatLon.value()->latLonInterp.value(),
                                                                                        resol_);
      size_t numstates = incrementsToInterp_.size();
      for (size_t jm=0; jm < numstates; jm++) {
        // This supports output on model levels only; to output on pressure levels would need to
        // read in a reference background from which to read the vertical pressure coordinate
        latLonWriter_.interpolateAndWrite(incrementsToInterp_[jm]);
        Log::test() << latLonWriter_ << std::endl;
      }
    }

    return 0;
  }
// -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    ConvertToLatLonParameters_ params;
    params.outputSchema(outputPath);
  }
// -----------------------------------------------------------------------------
  void validateConfig(const eckit::Configuration & fullConfig) const override {
    ConvertToLatLonParameters_ params;
    params.validate(fullConfig);
  }
// -------------------------------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::ConvertToLatLon<" + MODEL::name() + ">";
  }
};

}  // namespace oops

#endif  // OOPS_RUNS_CONVERTTOLATLON_H_
