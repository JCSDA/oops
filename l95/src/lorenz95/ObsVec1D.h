/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_OBSVEC1D_H_
#define LORENZ95_OBSVEC1D_H_

#include <ostream>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "lorenz95/ObsTable.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace lorenz95 {
  class GomL95;
  class ObsTable;

// -----------------------------------------------------------------------------
/// Vector in observation space
/*!
 *  ObsVec1D is implemented as an STL vector
 */

class ObsVec1D : public util::Printable,
                 private util::ObjectCounter<ObsVec1D> {
 public:
  static const std::string classname() {return "lorenz95::ObsVec1D";}

  explicit ObsVec1D(const ObsTable &);
  ObsVec1D(const ObsVec1D &, const bool copy = true);
  ~ObsVec1D() {}

  ObsVec1D & operator= (const ObsVec1D &);
  ObsVec1D & operator*= (const double &);
  ObsVec1D & operator+= (const ObsVec1D &);
  ObsVec1D & operator-= (const ObsVec1D &);
  ObsVec1D & operator*= (const ObsVec1D &);
  ObsVec1D & operator/= (const ObsVec1D &);

  void zero();
  void axpy(const double &, const ObsVec1D &);
  void invert();
  void random();
  double dot_product_with(const ObsVec1D &) const;
  double rms() const;

  unsigned int size() const {return data_.size();}
  double & operator() (const unsigned int ii) {return data_[ii];}
  const double & operator() (const unsigned int ii) const {return data_[ii];}

// I/O
  void read(const std::string &);
  void save(const std::string &) const;

 private:
  void print(std::ostream &) const;

  const ObsTable & obsdb_;
  std::vector<double> data_;
};
//-----------------------------------------------------------------------------
}  // namespace lorenz95
#endif  // LORENZ95_OBSVEC1D_H_
