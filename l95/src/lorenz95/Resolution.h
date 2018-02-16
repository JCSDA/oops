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

#include "eckit/config/Configuration.h"
#include "util/Printable.h"

namespace lorenz95 {

// -----------------------------------------------------------------------------
/// Handles resolution.

class Resolution : public util::Printable {
 public:
  explicit Resolution(const eckit::Configuration & conf): resol_(conf.getInt("resol")) {}
  explicit Resolution(const int resol): resol_(resol) {}
  Resolution(const Resolution & other): resol_(other.resol_) {}
  ~Resolution() {}

  std::vector<double> getLats() const;
  std::vector<double> getLons() const;
  std::vector<double> getLevs() const;
  std::vector<int> getMask(const int &) const;

  int npoints() const {return resol_;}

 private:
  Resolution & operator=(const Resolution &);
  void print(std::ostream & os) const {os << resol_;}
  const int resol_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_RESOLUTION_H_
