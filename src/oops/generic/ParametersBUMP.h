/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_GENERIC_PARAMETERSBUMP_H_
#define OOPS_GENERIC_PARAMETERSBUMP_H_

#include <sstream>
#include <string>
#include <vector>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/Ensemble.h"
#include "oops/base/EnsemblesCollection.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/Variables.h"
#include "oops/generic/oobump_f.h"
#include "oops/generic/UnstructuredGrid.h"
#include "oops/interface/ErrorCovariance.h"
#include "oops/interface/LocalizationBase.h"
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
  typedef ErrorCovariance<MODEL>                ErrorCovariance_;
  typedef Geometry<MODEL>                       Geometry_;
  typedef Increment<MODEL>                      Increment_;
  typedef State<MODEL>                          State_;

  typedef Ensemble<MODEL>                       Ensemble_;
  typedef boost::shared_ptr<Ensemble<MODEL> >   EnsemblePtr_;
  typedef EnsemblesCollection<MODEL>            EnsemblesCollection_;

 public:
  static const std::string classname() {return "oops::ParametersBUMP";}
  explicit ParametersBUMP(const eckit::Configuration &);
  ~ParametersBUMP();

  void estimate() const;
  void write() const;

 private:
  const eckit::LocalConfiguration conf_;
  int colocated_;
  int keyBUMP_;
};

// =============================================================================

template<typename MODEL>
ParametersBUMP<MODEL>::ParametersBUMP(const eckit::Configuration & conf)
  : conf_(conf), colocated_(1), keyBUMP_(0)
{
  Log::trace() << "ParametersBUMP<MODEL>::ParametersBUMP construction starting" << std::endl;
  util::Timer timer(classname(), "ParametersBUMP");

  const eckit::Configuration * fconf = &conf;

//  Setup resolution
  const eckit::LocalConfiguration resolConfig(conf_, "resolution");
  const Geometry_ resol(resolConfig);

// Setup variables
  const eckit::LocalConfiguration varConfig(conf_, "variables");
  const Variables vars(varConfig);

// Setup time
  const util::DateTime date(conf_.getString("date"));

// Setup background state
  const eckit::LocalConfiguration backgroundConfig(conf_, "background");
  State_ xx(resol, vars, backgroundConfig);

// Setup dummy increment
  Increment_ dx(resol, vars, date);

// Define unstructured grid coordinates
  colocated_ = 1;  // conf_.getString("colocated") TODO
  UnstructuredGrid ug;
  dx.ug_coord(ug, colocated_);

// Compute the ensemble of perturbations at time of xx
  const eckit::LocalConfiguration ensembleConfig(conf_, "ensemble");
  EnsemblePtr_ ens_ptr(new Ensemble_(xx.validTime(), ensembleConfig));
  ens_ptr->linearize(xx, xx, resol);
  EnsemblesCollection_::getInstance().put(xx.validTime(), ens_ptr);
  int ens1_ne = ens_ptr->size();

// Setup pseudo ensemble
  boost::ptr_vector<Increment_> bkgens;
  int ens2_ne = 0;
  if (conf_.has("covariance")) {
  // Setup background error covariance
    const eckit::LocalConfiguration covarConfig(conf_, "covariance");
    ErrorCovariance_ cov(resol, vars, covarConfig, xx, xx);

  // Compute a pseudo ensemble using randomization
    ens2_ne = covarConfig.getInt("pseudoens_size");
    for (int ii = 0; ii < ens2_ne; ++ii) {
      Increment_ * incr = new Increment_(resol, vars, date);
      cov.randomize(*incr);
      bkgens.push_back(incr);
    }
  }

// Create BUMP
  create_oobump_f90(keyBUMP_, ug.toFortran(), &fconf, ens1_ne, 1, ens2_ne, 1);

// Copy ensemble members
  const double rk = sqrt((static_cast<double>(ens1_ne) - 1.0));
  for (int ie = 0; ie < ens1_ne; ++ie) {
    Log::info() << "Copy ensemble member " << ie+1 << " / "
                 << ens1_ne << " to BUMP" << std::endl;

    // Renormalize member
    dx = (*ens_ptr)[ie];
    dx *= rk;

    // Define unstructured grid field
    UnstructuredGrid ugmem;
    dx.field_to_ug(ugmem, colocated_);

    // Copy field into BUMP ensemble
    add_oobump_member_f90(keyBUMP_, ugmem.toFortran(), ie+1, bump::readEnsMember);
  }

// Copy pseudo-ensemble members
  for (int ie = 0; ie < ens2_ne; ++ie) {
     Log::info() << "Copy pseudo ensemble member " << ie+1 << " / "
                 << ens2_ne << " to BUMP" << std::endl;

    // Define unstructured grid field
    UnstructuredGrid ugmem;
    bkgens[ie].field_to_ug(ugmem, colocated_);

    // Copy field into BUMP ensemble
    add_oobump_member_f90(keyBUMP_, ugmem.toFortran(), ie+1, bump::readPseudoEnsMember);
  }

  Log::trace() << "ParametersBUMP:ParametersBUMP constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ParametersBUMP<MODEL>::~ParametersBUMP() {
  Log::trace() << "ParametersBUMP<MODEL>::~ParametersBUMP destruction starting" << std::endl;
  util::Timer timer(classname(), "~ParametersBUMP");
  delete_oobump_f90(keyBUMP_);
  Log::trace() << "ParametersBUMP:~ParametersBUMP destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ParametersBUMP<MODEL>::estimate() const {
  Log::trace() << "ParametersBUMP::estimate starting" << std::endl;
  util::Timer timer(classname(), "estimate");

//  Setup resolution
  const eckit::LocalConfiguration resolConfig(conf_, "resolution");
  const Geometry_ resol(resolConfig);

// Setup variables
  const eckit::LocalConfiguration varConfig(conf_, "variables");
  const Variables vars(varConfig);

// Setup time
  const util::DateTime date(conf_.getString("date"));

// Setup dummy increment
  Increment_ dx(resol, vars, date);
  dx.zero();

// Setup unstructured grid
  UnstructuredGrid ug;

// Read data from files
  if (conf_.has("input")) {
    std::vector<eckit::LocalConfiguration> inputConfigs;
    conf_.get("input", inputConfigs);
    for (const auto & conf : inputConfigs) {
      dx.read(conf);
      dx.field_to_ug(ug, colocated_);
      std::string param = conf.getString("parameter");
      const int nstr = param.size();
      const char *cstr = param.c_str();
      set_oobump_param_f90(keyBUMP_, nstr, cstr, ug.toFortran());
    }
  }

// Estimate parameters
  run_oobump_drivers_f90(keyBUMP_);

  Log::trace() << "ParametersBUMP:estimate done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ParametersBUMP<MODEL>::write() const {
  Log::trace() << "ParametersBUMP::write starting" << std::endl;
  util::Timer timer(classname(), "write");

//  Setup resolution
  const eckit::LocalConfiguration resolConfig(conf_, "resolution");
  const Geometry_ resol(resolConfig);

// Setup variables
  const eckit::LocalConfiguration varConfig(conf_, "variables");
  const Variables vars(varConfig);

// Setup time
  const util::DateTime date(conf_.getString("date"));

// Setup dummy increment
  Increment_ dx(resol, vars, date);
  dx.zero();

// Setup unstructured grid
  UnstructuredGrid ug;
  dx.field_to_ug(ug, colocated_);

// Write parameters
  std::vector<eckit::LocalConfiguration> outputConfigs;
  conf_.get("output", outputConfigs);
  for (const auto & conf : outputConfigs) {
    std::string param = conf.getString("parameter");
    const int nstr = param.size();
    const char *cstr = param.c_str();
    get_oobump_param_f90(keyBUMP_, nstr, cstr, ug.toFortran());
    dx.field_from_ug(ug);
    dx.write(conf);
    Log::test() << "Norm of " << param << ": " << dx.norm() << std::endl;
  }
  Log::trace() << "ParametersBUMP::write done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_PARAMETERSBUMP_H_
