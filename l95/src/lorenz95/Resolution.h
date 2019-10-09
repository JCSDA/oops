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
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"
#include "lorenz95/Iterator.h"
#include "oops/util/Printable.h"

namespace lorenz95 {

class Iterator;

// -----------------------------------------------------------------------------
/// Handles resolution.

class Resolution : public util::Printable {
 public:
  Resolution(const eckit::Configuration & conf, const eckit::mpi::Comm & comm)
              : resol_(conf.getInt("resol")), comm_(comm) {}
  explicit Resolution(const int resol): resol_(resol), comm_(eckit::mpi::comm("self")) {}
  Resolution(const Resolution & other): resol_(other.resol_), comm_(other.comm_) {}
  ~Resolution() {}

  int npoints() const {return resol_;}

  Iterator begin() const;
  Iterator end() const;
  const eckit::mpi::Comm & getComm() const {return comm_;}

 private:
  Resolution & operator=(const Resolution &);
  void print(std::ostream & os) const {os << resol_;}
  const int resol_;
  const eckit::mpi::Comm & comm_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_RESOLUTION_H_
