/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_GENERIC_LOCALIZATIONNICAS_H_
#define OOPS_GENERIC_LOCALIZATIONNICAS_H_

#include <sstream>
#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/generic/nicas_f.h"
#include "oops/generic/UnstructuredGrid.h"
#include "oops/interface/LocalizationBase.h"
#include "util/Logger.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// NICAS localization matrix.

template<typename MODEL>
class LocalizationNICAS : public LocalizationBase<MODEL> {
  typedef typename MODEL::Increment             Increment_;
  typedef typename MODEL::State                 State_;

 public:
  LocalizationNICAS(const State_ &, const eckit::Configuration &);
  ~LocalizationNICAS();

  void multiply(Increment_ &) const;

 private:
  void print(std::ostream &) const;

  int keyNicas_;
};

// =============================================================================

template<typename MODEL>
LocalizationNICAS<MODEL>::LocalizationNICAS(const State_ & xx, const eckit::Configuration & conf) {
  const eckit::Configuration * fconf = &conf;

// Get lat/lon/mask from the unstructured grid
  UnstructuredGrid ugrid;
  xx.convert_to(ugrid);
  std::vector<double> lats = ugrid.getLats();
  std::vector<double> lons = ugrid.getLons();
  std::vector<double> levs = ugrid.getLevs();
  std::vector<int> cmask;
  for (int jlev = 0; jlev < levs.size(); ++jlev) {
    std::vector<int> tmp = ugrid.getCmask(jlev);
    cmask.insert(cmask.end(), tmp.begin(), tmp.end());
  }
  int nh = lats.size();
  int nv = levs.size();
  create_nicas_f90(keyNicas_, &fconf, nh, &lats[0], &lons[0], nv, &levs[0], &cmask[0]);
  Log::trace() << "LocalizationNICAS:LocalizationNICAS constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
LocalizationNICAS<MODEL>::~LocalizationNICAS() {
  delete_nicas_f90(keyNicas_);
  Log::trace() << "LocalizationNICAS:~LocalizationNICAS destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationNICAS<MODEL>::multiply(Increment_ & dx) const {
  Log::trace() << "LocalizationNICAS:multiply starting" << std::endl;

  UnstructuredGrid ugrid;
  dx.convert_to(ugrid);

  nicas_multiply_f90(keyNicas_, ugrid.toFortran());

  dx.convert_from(ugrid);

  Log::trace() << "LocalizationNICAS:multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationNICAS<MODEL>::print(std::ostream & os) const {
  os << "LocalizationNICAS<MODEL>::print not implemeted yet";
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_LOCALIZATIONNICAS_H_
