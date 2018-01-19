/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_GENERIC_LOCALIZATIONHDIAG_NICAS_H_
#define OOPS_GENERIC_LOCALIZATIONHDIAG_NICAS_H_

#include <sstream>
#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/Ensemble.h"
#include "oops/base/EnsemblesCollection.h"
#include "oops/base/Variables.h"
#include "oops/generic/hdiag_nicas_f.h"
#include "oops/generic/UnstructuredGrid.h"
#include "oops/interface/LocalizationBase.h"
#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/Logger.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// HDIAG_NICAS localization matrix.

template<typename MODEL>
class LocalizationHDIAG_NICAS : public LocalizationBase<MODEL> {
  typedef typename MODEL::Geometry              Geometry_;
  typedef typename MODEL::Increment             Increment_;
  typedef typename MODEL::State                 State_;

  typedef boost::shared_ptr<Ensemble<MODEL> >   EnsemblePtr_;
  typedef EnsemblesCollection<MODEL>            EnsemblesCollection_;

 public:
  LocalizationHDIAG_NICAS(const Geometry_ &, const eckit::Configuration &);
  ~LocalizationHDIAG_NICAS();

  void multiply(Increment_ &) const;

 private:
  void print(std::ostream &) const;

  int keyHdiag_nicas_;
};

// =============================================================================

template<typename MODEL>
LocalizationHDIAG_NICAS<MODEL>::LocalizationHDIAG_NICAS(const Geometry_ & resol, const eckit::Configuration & conf) {
  const eckit::Configuration * fconf = &conf;

// Setup variables
  const eckit::LocalConfiguration varConfig(conf, "variables");
  const Variables vars(varConfig);

// Setup dummy time
  const util::DateTime date;

// Setup dummy increment
  Increment_ dx(resol, vars, date);

// Convert to unstructured grid
  UnstructuredGrid ug;
//  dx.define(ug);

// Get sizes and coordinates from the unstructured grid
  int nc0a = ug.getSize(1);
  int nl0 = ug.getSize(2);
  int nv = ug.getSize(3);
  int nts = 1; //ug.getSize(4); not read yet for 4D
  std::vector<double> lon = ug.getLon();
  std::vector<double> lat = ug.getLat();
  std::vector<double> area = ug.getArea();
  std::vector<double> vunit = ug.getVunit();
  std::vector<int> imask = ug.getImask();

  int ens1_ne;
  std::vector<double> ens1;
  int new_hdiag = conf.getInt("new_hdiag");
  if (new_hdiag==1) {
    std::vector<eckit::LocalConfiguration> confs;
    conf.get("hdiag_ensemble", confs);
    ens1_ne = confs.size();
    Log::info() << "HDIAG ensemble: " << ens1_ne << " members" << std::endl;
    ASSERT(confs.size()==ens1_ne*nts);
    for (unsigned int ie = 0; ie < ens1_ne; ++ie) {
      Log::info() << "Read member " << ie+1 << std::endl;
      State_ xmem(resol,confs[ie]);
      Log::info() << "Convert member to unstructured grid" << std::endl;
      xmem.convert_to(ug);
      std::vector<double> tmp = ug.getData();
      ens1.insert(ens1.end(), tmp.begin(), tmp.end());
    }
  }
  else
  {
    ens1_ne = 4;
    for (unsigned int ie = 0; ie < ens1_ne; ++ie) {
      std::vector<double> tmp(nc0a*nl0*nv*nts);
      std::fill(tmp.begin(),tmp.end(),-999.0);
      ens1.insert(ens1.end(), tmp.begin(), tmp.end());
    }
  }
  create_hdiag_nicas_f90(keyHdiag_nicas_, &fconf, nc0a, nl0, nv, nts, ens1_ne, &lon[0], &lat[0], &area[0], &vunit[0], &imask[0], &ens1[0]);
  Log::trace() << "LocalizationHDIAG_NICAS:LocalizationHDIAG_NICAS constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
LocalizationHDIAG_NICAS<MODEL>::~LocalizationHDIAG_NICAS() {
  delete_hdiag_nicas_f90(keyHdiag_nicas_);
  Log::trace() << "LocalizationHDIAG_NICAS:~LocalizationHDIAG_NICAS destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationHDIAG_NICAS<MODEL>::multiply(Increment_ & dx) const {
  Log::trace() << "LocalizationHDIAG_NICAS:multiply starting" << std::endl;

  UnstructuredGrid ug;

  MPI_Barrier(MPI_COMM_WORLD);
  boost::posix_time::ptime ti0 = boost::posix_time::microsec_clock::local_time();
//  dx.define(ug);
  boost::posix_time::ptime t0 = boost::posix_time::microsec_clock::local_time();
  boost::posix_time::time_duration diff0 = t0 - ti0;
  Log::info() << "define time: " << diff0.total_nanoseconds()/1000 << " microseconds" << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);
  boost::posix_time::ptime ti1 = boost::posix_time::microsec_clock::local_time();
  dx.convert_to(ug);
  boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
  boost::posix_time::time_duration diff1 = t1 - ti1;
  Log::info() << "convert_to time: " << diff1.total_nanoseconds()/1000 << " microseconds" << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);
  boost::posix_time::ptime ti2 = boost::posix_time::microsec_clock::local_time();
  hdiag_nicas_multiply_f90(keyHdiag_nicas_, ug.toFortran());
  boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
  boost::posix_time::time_duration diff2 = t2 - ti2;
  Log::info() << "multiply time: " << diff2.total_nanoseconds()/1000 << " microseconds" << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);
  boost::posix_time::ptime ti3 = boost::posix_time::microsec_clock::local_time();
  dx.convert_from(ug);
  boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
  boost::posix_time::time_duration diff3 = t3 - ti3;
  Log::info() << "convert_from time: " << diff3.total_nanoseconds()/1000 << " microseconds" << std::endl;

  Log::trace() << "LocalizationHDIAG_NICAS:multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationHDIAG_NICAS<MODEL>::print(std::ostream & os) const {
  os << "LocalizationHDIAG_NICAS<MODEL>::print not implemeted yet";
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_LOCALIZATIONHDIAG_NICAS_H_
