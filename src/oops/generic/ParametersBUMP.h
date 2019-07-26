/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_PARAMETERSBUMP_H_
#define OOPS_GENERIC_PARAMETERSBUMP_H_

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/date_time/posix_time/posix_time.hpp>

#include "eckit/config/Configuration.h"

#include "oops/assimilation/Increment4D.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/Variables.h"
#include "oops/generic/oobump_f.h"
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
  typedef Geometry<MODEL>                         Geometry_;
  typedef Increment<MODEL>                        Increment_;
  typedef Increment4D<MODEL>                      Increment4D_;
  typedef State<MODEL>                            State_;
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

  int get_bump() const {return keyBUMP_;}
  void clean_bump() const {delete_oobump_f90(keyBUMP_);}
  void write() const;

 private:
  const Geometry_ resol_;
  const Variables vars_;
  std::vector<util::DateTime> timeslots_;
  const eckit::LocalConfiguration conf_;
  int keyBUMP_;
};

// =============================================================================

template<typename MODEL>
ParametersBUMP<MODEL>::ParametersBUMP(const Geometry_ & resol,
                                      const Variables & vars,
                                      const std::vector<util::DateTime> & timeslots,
                                      const eckit::Configuration & conf,
                                      const EnsemblePtr_ ens,
                                      const EnsemblePtr_ pseudo_ens)
  : resol_(resol), vars_(vars), timeslots_(timeslots), conf_(conf), keyBUMP_(0)
{
  Log::trace() << "ParametersBUMP<MODEL>::ParametersBUMP construction starting" << std::endl;
  util::Timer timer(classname(), "ParametersBUMP");

// Setup BUMP configuration
  const eckit::LocalConfiguration BUMPConfig(conf_, "bump");
  const eckit::Configuration * fconf = &BUMPConfig;

// Setup colocation
  int colocated = 1;
  if (BUMPConfig.has("colocated")) colocated = BUMPConfig.getInt("colocated");

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
  create_oobump_f90(keyBUMP_, ug.toFortran(), &fconf, ens1_ne, 1, ens2_ne, 1);

// Add ensemble members
  Log::info() << "Add ensemble members" << std::endl;
  for (int ie = 0; ie < ens1_ne; ++ie) {
    Log::info() << "   Copy ensemble member " << ie+1 << " / "
                << ens1_ne << " to BUMP" << std::endl;

  // Copy member
    dx = (*ens)[ie];

  // Renormalize member
    const double rk = sqrt((static_cast<double>(ens1_ne) - 1.0));
    dx *= rk;

  // Define unstructured grid field
    dx.field_to_ug(ug);

  // Copy field into BUMP ensemble
    add_oobump_member_f90(keyBUMP_, ug.toFortran(), ie+1, bump::readEnsMember);
  }

// Add pseudo-ensemble members
  Log::info() << "Add pseudo-ensemble members" << std::endl;
  for (int ie = 0; ie < ens2_ne; ++ie) {
    Log::info() << "   Copy pseudo-ensemble member " << ie+1 << " / "
                << ens2_ne << " to BUMP" << std::endl;

  // Copy member
    dx = (*pseudo_ens)[ie];

  // Define unstructured grid field
    dx.field_to_ug(ug);

  // Copy field into BUMP pseudo-ensemble
    add_oobump_member_f90(keyBUMP_, ug.toFortran(), ie+1, bump::readPseudoEnsMember);
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
      const int nstr = param.size();
      const char *cstr = param.c_str();
      set_oobump_param_f90(keyBUMP_, nstr, cstr, ug.toFortran());
    }
  }

// Estimate parameters
  run_oobump_drivers_f90(keyBUMP_);

// Copy BUMP test file
  const std::string bump_test = BUMPConfig.getString("prefix") + ".test.0000";
  std::ifstream infile(bump_test);
  std::string line;
//  while (std::getline(infile, line)) Log::test() << line << std::endl;

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
  int colocated;
  get_oobump_colocated_f90(keyBUMP_, colocated);
  UnstructuredGrid ug(colocated, timeslots_.size());
  dx.ug_coord(ug);

// Write parameters
  std::vector<eckit::LocalConfiguration> outputConfigs;
  conf_.get("output", outputConfigs);
  for (const auto & conf : outputConfigs) {
  // Get parameter from BUMP
    std::string param = conf.getString("parameter");
    const int nstr = param.size();
    const char *cstr = param.c_str();
    get_oobump_param_f90(keyBUMP_, nstr, cstr, ug.toFortran());
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
