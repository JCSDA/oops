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
  typedef typename MODEL::Geometry              Geometry_;
  typedef typename MODEL::Increment             Increment_;

 public:
  LocalizationNICAS(const Geometry_ &, const eckit::Configuration &);
  ~LocalizationNICAS();

  void multiply(Increment_ &) const;

 private:
  void print(std::ostream &) const;

  int keyNicas_;
};

// =============================================================================

template<typename MODEL>
LocalizationNICAS<MODEL>::LocalizationNICAS(const Geometry_ & grid, const eckit::Configuration & conf) {
  const eckit::Configuration * fconf = &conf;
  std::vector<int> dims = grid.getDims();
  std::vector<double> lats = grid.getLats();
  std::vector<double> lons = grid.getLons();
  std::vector<double> levs = grid.getLevs();
  std::vector<double> area = grid.getArea();
  std::vector<int> mask;
  for (int jlev = 0; jlev < levs.size(); ++jlev) {
    std::vector<int> tmp = grid.getMask(jlev);
    mask.insert(mask.end(), tmp.begin(), tmp.end());
  }
  int ndims = dims.size();
  int nh = lats.size();
  int nv = levs.size();
  create_nicas_f90(keyNicas_, &fconf, ndims, &dims[0], nh, &lats[0], &lons[0], nv, &levs[0], &area[0], &mask[0]);
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
