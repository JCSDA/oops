/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_PARAMETERSBUMP_H_
#define OOPS_GENERIC_PARAMETERSBUMP_H_

#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <boost/date_time/posix_time/posix_time.hpp>

#include "eckit/config/Configuration.h"

#include "oops/assimilation/Increment4D.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/Variables.h"
#include "oops/generic/OoBump.h"
#include "oops/generic/UnstructuredGrid.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// BUMP diagnostics.

template<typename MODEL>
class ParametersBUMP {
  typedef Geometry<MODEL>                             Geometry_;
  typedef Increment<MODEL>                            Increment_;
  typedef Increment4D<MODEL>                          Increment4D_;
  typedef State<MODEL>                                State_;
  typedef boost::shared_ptr<IncrementEnsemble<MODEL>> EnsemblePtr_;

 public:
  static const std::string classname() {return "oops::ParametersBUMP";}
  ParametersBUMP(const Geometry_ &,
                 const Variables &,
                 const std::vector<util::DateTime> &,
                 const eckit::Configuration &,
                 const EnsemblePtr_ ens = NULL,
                 const EnsemblePtr_ pseudo_ens = NULL);
  ~ParametersBUMP();

  OoBump & getOoBump() {return *ooBump_;}
  void write() const;

 private:
  const Geometry_ resol_;
  const Variables vars_;
  std::vector<util::DateTime> timeslots_;
  const eckit::LocalConfiguration conf_;
  std::unique_ptr<OoBump> ooBump_;
};

// =============================================================================

template<typename MODEL>
ParametersBUMP<MODEL>::ParametersBUMP(const Geometry_ & resol,
                                      const Variables & vars,
                                      const std::vector<util::DateTime> & timeslots,
                                      const eckit::Configuration & conf,
                                      const EnsemblePtr_ ens,
                                      const EnsemblePtr_ pseudo_ens)
  : resol_(resol), vars_(vars), timeslots_(timeslots), conf_(conf), ooBump_()
{
  Log::trace() << "ParametersBUMP<MODEL>::ParametersBUMP construction starting" << std::endl;
  util::Timer timer(classname(), "ParametersBUMP");

// Setup BUMP configuration
  const eckit::LocalConfiguration BUMPConfig(conf_, "bump");

// Setup colocation
  int colocated = 1;
  if (BUMPConfig.has("colocated")) colocated = BUMPConfig.getInt("colocated");

// Setup members release
  int release_members = 0;
  if (BUMPConfig.has("release_members")) release_members = BUMPConfig.getInt("release_members");

// Setup dummy increment
  Increment4D_ dx(resol_, vars_, timeslots_);

// Define unstructured grid coordinates
  UnstructuredGrid ug(colocated, timeslots_.size());
  dx.ug_coord(ug);

// Get ensemble size if ensemble is available
  int ens1_ne = 0;
  if (ens) ens1_ne = ens->size();

// Get pseudo-ensemble size if pseudo-ensemble is available
  int ens2_ne = 0;
  if (pseudo_ens) ens2_ne = pseudo_ens->size();

// Create BUMP
  Log::info() << "Create BUMP" << std::endl;
  ooBump_.reset(new OoBump(ug, BUMPConfig, ens1_ne, 1, ens2_ne, 1));

// Transfer/copy ensemble members to BUMP
  if (release_members == 1) {
    Log::info() << "Transfer ensemble members to BUMP" << std::endl;
  } else {
    Log::info() << "Copy ensemble members to BUMP" << std::endl;
  }
  for (int ie = 0; ie < ens1_ne; ++ie) {
    Log::info() << "   Member " << ie+1 << " / " << ens1_ne << std::endl;;

  // Copy member
    if (release_members == 1) {
      dx = (*ens)[0];
    } else {
      dx = (*ens)[ie];
    }

  // Define unstructured grid
    dx.field_to_ug(ug);

  // Copy data to BUMP
    ooBump_->addMember(ug, ie);

    if (release_members == 1) {
    // Release ensemble member
      ens->releaseMember();
    }
  }

// Transfer/copy pseudo-ensemble members to BUMP
  if (release_members == 1) {
    Log::info() << "Transfer pseudo-ensemble members to BUMP" << std::endl;
  } else {
    Log::info() << "Copy pseudo-ensemble members to BUMP" << std::endl;
  }
  for (int ie = 0; ie < ens2_ne; ++ie) {
    Log::info() << "   Member " << ie+1 << " / " << ens2_ne << std::endl;

  // Copy member
    if (release_members == 1) {
      dx = (*pseudo_ens)[0];
    } else {
      dx = (*pseudo_ens)[ie];
    }

  // Define unstructured grid
    dx.field_to_ug(ug);

  // Copy data to BUMP
    ooBump_->addPseudoMember(ug, ie);

    if (release_members == 1) {
    // Release ensemble member
      pseudo_ens->releaseMember();
    }
  }

// Read data from files
  Log::info() << "Read data from files" << std::endl;
  if (conf_.has("input")) {
  // Set BUMP input parameters
    std::vector<eckit::LocalConfiguration> inputConfigs;
    conf_.get("input", inputConfigs);
    for (const auto & conf : inputConfigs) {
    // Read parameter for the specified timeslot
      const util::DateTime date(conf.getString("date"));
      bool found = false;
      dx.zero();
      for (unsigned jsub = 0; jsub < timeslots_.size(); ++jsub) {
        if (date == timeslots_[jsub]) {
          found = true;
          dx[dx.first()+jsub].read(conf);
        }
      }
      ASSERT(found);
      dx.field_to_ug(ug);

    // Set parameter to BUMP
      std::string param = conf.getString("parameter");
      ooBump_->setParam(param, ug);
    }
  }

// Estimate parameters
  ooBump_->runDrivers();

  if (release_members == 1) {
  // Transfer ensemble members from BUMP
    Log::info() << "Transfer ensemble members from BUMP" << std::endl;
    for (int ie = 0; ie < ens1_ne; ++ie) {
      Log::info() << "   Member " << ie+1 << " / " << ens1_ne << std::endl;

    // Copy data from BUMP
      ooBump_->removeMember(ug, ie);

    // Reset ensemble member
      dx.field_from_ug(ug);
      ens->resetMember(dx);
    }

  // Transfer pseudo-ensemble members from BUMP
    Log::info() << "Transfer pseudo-ensemble members from BUMP" << std::endl;
    for (int ie = 0; ie < ens2_ne; ++ie) {
      Log::info() << "   Member " << ie+1 << " / " << ens2_ne << std::endl;

    // Copy data from BUMP
      ooBump_->removePseudoMember(ug, ie);

    // Reset pseudo-ensemble member
      dx.field_from_ug(ug);
      pseudo_ens->resetMember(dx);
    }
  }

  Log::trace() << "ParametersBUMP:ParametersBUMP constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ParametersBUMP<MODEL>::~ParametersBUMP() {
  Log::trace() << "ParametersBUMP<MODEL>::~ParametersBUMP destruction starting" << std::endl;
  util::Timer timer(classname(), "~ParametersBUMP");
  Log::trace() << "ParametersBUMP:~ParametersBUMP destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ParametersBUMP<MODEL>::write() const {
  Log::trace() << "ParametersBUMP::write starting" << std::endl;
  util::Timer timer(classname(), "write");


// Setup dummy increment
  Increment4D_ dx(resol_, vars_, timeslots_);
  dx.zero();

// Setup unstructured grid
  UnstructuredGrid ug(ooBump_->getColocated(), timeslots_.size());
  dx.ug_coord(ug);

// Write parameters
  std::vector<eckit::LocalConfiguration> outputConfigs;
  conf_.get("output", outputConfigs);
  for (const auto & conf : outputConfigs) {
  // Get parameter from BUMP
    std::string param = conf.getString("parameter");
    ooBump_->getParam(param, ug);
    dx.field_from_ug(ug);

  // Write parameter for the specified timeslot
    const util::DateTime date(conf.getString("date"));
    bool found = false;
    for (unsigned jsub = 0; jsub < timeslots_.size(); ++jsub) {
      int isub = jsub+dx.first();
      if (date == timeslots_[jsub]) {
        found = true;
        dx[isub].write(conf);
        Log::test() << "Norm of " << param << " at " << date << ": " << std::scientific
                    << std::setprecision(3) << dx[isub].norm() << std::endl;
      }
    }
    ASSERT(found);
  }
  Log::trace() << "ParametersBUMP::write done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_PARAMETERSBUMP_H_
