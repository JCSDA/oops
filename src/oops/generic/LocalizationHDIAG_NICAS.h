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

// Get lat/lon/area/vunit/mask from the unstructured grid
  UnstructuredGrid ugrid;
  dx.convert_to(ugrid); // TODO: add key, not to copy data
  std::vector<double> lats = ugrid.getLats();
  std::vector<double> lons = ugrid.getLons();
  std::vector<double> areas = ugrid.getAreas();
  std::vector<double> vunit = ugrid.getVunit();
  std::vector<int> mask3d;
  for (int jlev = 0; jlev < vunit.size(); ++jlev) {
    std::vector<int> tmp = ugrid.getMask3d(jlev);
    mask3d.insert(mask3d.end(), tmp.begin(), tmp.end());
  }
  std::vector<int> mask2d = ugrid.getMask2d();
  std::vector<int> glbind = ugrid.getGlbInd();
  int nv = ugrid.getNvar3d();
  int nc0a = lats.size();
  int nlev = vunit.size();

  int nts = 1;
  std::vector<eckit::LocalConfiguration> confs;
  conf.get("hdiag_ensemble", confs);
  int ens1_ne = confs.size();
  Log::info() << "HDIAG ensemble: "  << ens1_ne << " members / " << nts << " timeslots" << std::endl;
  ASSERT(confs.size()==ens1_ne*nts);
  std::vector<double> ens1;
  int i=0;
  for (unsigned int ie = 0; ie < ens1_ne; ++ie) {
    for (unsigned int its = 0; its < nts; ++its) { 
      Log::info() << "Read member " << ie+1 << " at timeslot " << its+1 << " (total index " << i << ")" << std::endl;
      State_ xmem(resol,confs[i]);
      Log::info() << "Convert member to unstructured grid" << std::endl;
      UnstructuredGrid umem;
      xmem.convert_to(umem);
      std::vector<double> tmp = umem.getData();
      ens1.insert(ens1.end(), tmp.begin(), tmp.end());
      i++;
    }
  }
  create_hdiag_nicas_f90(keyHdiag_nicas_, &fconf, nv, nc0a, &lats[0], &lons[0], &areas[0], nlev, &vunit[0], &mask3d[0], &mask2d[0], &glbind[0], nts, ens1_ne, &ens1[0]);
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

  UnstructuredGrid ugrid;
  dx.convert_to(ugrid);

  hdiag_nicas_multiply_f90(keyHdiag_nicas_, ugrid.toFortran());

  dx.convert_from(ugrid);

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
