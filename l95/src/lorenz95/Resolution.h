/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_RESOLUTION_H_
#define LORENZ95_RESOLUTION_H_

#include <iostream>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/config/Configuration.h"

#include "lorenz95/Iterator.h"

#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/util/Printable.h"

namespace lorenz95 {

class Iterator;

// -----------------------------------------------------------------------------
/// Handles resolution.

class Resolution : public util::Printable {
 public:
  Resolution(const eckit::Configuration & conf, const eckit::mpi::Comm & comm)
    : resol_(conf.getInt("resol")), comm_(comm)
    {ASSERT(comm_.size() == 1);}
  explicit Resolution(const int resol) : resol_(resol), comm_(oops::mpi::myself())
    {ASSERT(comm_.size() == 1);}

  int npoints() const {return resol_;}

  Iterator begin() const;
  Iterator end() const;
  std::vector<double> verticalCoord(std::string &) const;
  std::vector<size_t> variableSizes(const oops::Variables &) const;
  bool levelsAreTopDown() const {return true;}
  const eckit::mpi::Comm & getComm() const {return comm_;}
  void latlon(std::vector<double> &, std::vector<double> &, const bool) const;
  const atlas::FunctionSpace & functionSpace() const {return noFunctionSpace_;}
  const atlas::FieldSet & fields() const {return noFields_;}
  int closestTask(const double, const double) const { return 0; }

 private:
  void print(std::ostream & os) const {os << resol_;}
  const int resol_;
  const eckit::mpi::Comm & comm_;
  atlas::FunctionSpace noFunctionSpace_;
  atlas::FieldSet noFields_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_RESOLUTION_H_
