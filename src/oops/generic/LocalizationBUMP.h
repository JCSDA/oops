/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_LOCALIZATIONBUMP_H_
#define OOPS_GENERIC_LOCALIZATIONBUMP_H_

#include <sstream>
#include <string>
#include <vector>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/Ensemble.h"
#include "oops/base/EnsemblesCollection.h"
#include "oops/base/Variables.h"
#include "oops/generic/bump_f.h"
#include "oops/generic/UnstructuredGrid.h"
#include "oops/interface/LocalizationBase.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// BUMP localization matrix.

template<typename MODEL>
class LocalizationBUMP : public LocalizationBase<MODEL> {
  typedef typename MODEL::Geometry              Geometry_;
  typedef typename MODEL::Increment             Increment_;
  typedef typename MODEL::State                 State_;

  typedef boost::shared_ptr<Ensemble<MODEL> >   EnsemblePtr_;
  typedef EnsemblesCollection<MODEL>            EnsemblesCollection_;

 public:
  LocalizationBUMP(const Geometry_ &, const eckit::Configuration &);
  ~LocalizationBUMP();

  void multiply(Increment_ &) const;

 private:
  void print(std::ostream &) const;

  int keyBUMP_;
};

// =============================================================================

template<typename MODEL>
LocalizationBUMP<MODEL>::LocalizationBUMP(const Geometry_ & resol,
                                          const eckit::Configuration & conf) {
  const eckit::Configuration * fconf = &conf;

// Setup variables
  const eckit::LocalConfiguration varConfig(conf, "variables");
  const Variables vars(varConfig);

// Setup dummy time
  const util::DateTime date(conf.getString("date"));

// Setup dummy increment
  Increment_ dx(resol, vars, date);

// Define unstructured grid coordinates
  UnstructuredGrid ug;
  dx.ug_coord(ug);

// Define ensemble size
  int new_hdiag = conf.getInt("new_hdiag");
  int ens1_ne;
  if (new_hdiag == 1) {
    EnsemblePtr_ ens_ptr = EnsemblesCollection_::getInstance()[dx.validTime()];
    ens1_ne = ens_ptr->size();
  } else {
    ens1_ne = 0;
  }

// Create BUMP
  create_bump_f90(keyBUMP_, ug.toFortran(), &fconf, ens1_ne);

// Copy ensemble members
  if (new_hdiag == 1) {
    EnsemblePtr_ ens_ptr = EnsemblesCollection_::getInstance()[dx.validTime()];
    for (int ie = 0; ie < ens1_ne; ++ie) {
       Log::info() << "Copy ensemble member " << ie+1 << " / "
                   << ens1_ne << " to BUMP" << std::endl;

      // Define unstructured grid field
      UnstructuredGrid ugmem;
      (*ens_ptr)[ie].field_to_ug(ugmem);

      // Copy field into BUMP ensemble
      add_bump_member_f90(keyBUMP_, ugmem.toFortran(), ie+1);
    }
  }

// Run BUMP
  run_bump_drivers_f90(keyBUMP_);

  Log::trace() << "LocalizationBUMP:LocalizationBUMP constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
LocalizationBUMP<MODEL>::~LocalizationBUMP() {
  delete_bump_f90(keyBUMP_);
  Log::trace() << "LocalizationBUMP:~LocalizationBUMP destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationBUMP<MODEL>::multiply(Increment_ & dx) const {
  Log::trace() << "LocalizationBUMP:multiply starting" << std::endl;

  UnstructuredGrid ug;

  MPI_Barrier(MPI_COMM_WORLD);
  boost::posix_time::ptime ti1 = boost::posix_time::microsec_clock::local_time();
  dx.field_to_ug(ug);
  boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
  boost::posix_time::time_duration diff1 = t1 - ti1;
  Log::info() << "field_to_ug time: " << diff1.total_nanoseconds()/1000
              << " microseconds" << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);
  boost::posix_time::ptime ti2 = boost::posix_time::microsec_clock::local_time();
  bump_multiply_f90(keyBUMP_, ug.toFortran());
  boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
  boost::posix_time::time_duration diff2 = t2 - ti2;
  Log::info() << "multiply time: " << diff2.total_nanoseconds()/1000
              << " microseconds" << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);
  boost::posix_time::ptime ti3 = boost::posix_time::microsec_clock::local_time();
  dx.field_from_ug(ug);
  boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
  boost::posix_time::time_duration diff3 = t3 - ti3;
  Log::info() << "field_from_ug time: " << diff3.total_nanoseconds()/1000
              << " microseconds" << std::endl;

  Log::trace() << "LocalizationBUMP:multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationBUMP<MODEL>::print(std::ostream & os) const {
  os << "LocalizationBUMP<MODEL>::print not implemeted yet";
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_LOCALIZATIONBUMP_H_
