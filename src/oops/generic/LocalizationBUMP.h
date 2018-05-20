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

  int keybump_;
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
  const util::DateTime date;

// Setup dummy increment
  Increment_ dx(resol, vars, date);

// Unstructured grid
  UnstructuredGrid ug;
  dx.convert_to(ug);

// Get sizes and coordinates from the unstructured grid
  int nmga = ug.getSize(1);
  int nl0 = ug.getSize(2);
  int nv = ug.getSize(3);
  int nts = 1;  // ug.getSize(4); not read yet for 4D
  std::vector<double> lon = ug.getLon();
  std::vector<double> lat = ug.getLat();
  std::vector<double> area = ug.getArea();
  std::vector<double> vunit = ug.getVunit();
  std::vector<int> imask = ug.getImask();

  int ens1_ne;
  std::vector<double> ens1;
  int new_hdiag = conf.getInt("new_hdiag");
  if (new_hdiag == 1) {
    std::vector<eckit::LocalConfiguration> confs;
    conf.get("bump_ensemble", confs);
    ens1_ne = confs.size();
    Log::info() << "BUMP ensemble: " << ens1_ne << " members" << std::endl;
    for (int ie = 0; ie < ens1_ne; ++ie) {
      Log::info() << "Read member " << ie+1 << std::endl;
      State_ xmem(resol, confs[ie]);
      Log::info() << "Convert member to unstructured grid" << std::endl;
      UnstructuredGrid ugmem;
      xmem.convert_to(ugmem);
      std::vector<double> tmp = ugmem.getData();
      ens1.insert(ens1.end(), tmp.begin(), tmp.end());
    }
  } else {
    ens1_ne = 4;
    for (int ie = 0; ie < ens1_ne; ++ie) {
      std::vector<double> tmp(nmga*nl0*nv*nts);
      std::fill(tmp.begin(), tmp.end(), -999.0);
      ens1.insert(ens1.end(), tmp.begin(), tmp.end());
    }
  }
  create_bump_f90(keybump_, &fconf, nmga, nl0, nv, nts, &lon[0], &lat[0], &area[0],
                  &vunit[0], &imask[0], ens1_ne, &ens1[0]);
  Log::trace() << "LocalizationBUMP:LocalizationBUMP constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
LocalizationBUMP<MODEL>::~LocalizationBUMP() {
  delete_bump_f90(keybump_);
  Log::trace() << "LocalizationBUMP:~LocalizationBUMP destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationBUMP<MODEL>::multiply(Increment_ & dx) const {
  Log::trace() << "LocalizationBUMP:multiply starting" << std::endl;

  UnstructuredGrid ug;

  MPI_Barrier(MPI_COMM_WORLD);
  boost::posix_time::ptime ti1 = boost::posix_time::microsec_clock::local_time();
  dx.convert_to(ug);
  boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
  boost::posix_time::time_duration diff1 = t1 - ti1;
  Log::info() << "convert_to time: " << diff1.total_nanoseconds()/1000
              << " microseconds" << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);
  boost::posix_time::ptime ti2 = boost::posix_time::microsec_clock::local_time();
  bump_multiply_f90(keybump_, ug.toFortran());
  boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
  boost::posix_time::time_duration diff2 = t2 - ti2;
  Log::info() << "multiply time: " << diff2.total_nanoseconds()/1000
              << " microseconds" << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);
  boost::posix_time::ptime ti3 = boost::posix_time::microsec_clock::local_time();
  dx.convert_from(ug);
  boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
  boost::posix_time::time_duration diff3 = t3 - ti3;
  Log::info() << "convert_from time: " << diff3.total_nanoseconds()/1000
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
