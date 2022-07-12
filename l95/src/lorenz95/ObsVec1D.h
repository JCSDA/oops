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

#include <Eigen/Dense>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace lorenz95 {
  class ObsTable;
  template <typename DATATYPE> class ObsData1D;

// -----------------------------------------------------------------------------
/// Vector in observation space
/*!
 *  ObsVec1D is implemented as an STL vector
 */

class ObsVec1D : public util::Printable,
                 private util::ObjectCounter<ObsVec1D> {
 public:
  static const std::string classname() {return "lorenz95::ObsVec1D";}

  explicit ObsVec1D(const ObsTable &, const std::string & name = "");
  ObsVec1D(const ObsVec1D &);
  ~ObsVec1D() = default;

  ObsVec1D & operator= (const ObsVec1D &);
  ObsVec1D & operator*= (const double &);
  ObsVec1D & operator+= (const ObsVec1D &);
  ObsVec1D & operator-= (const ObsVec1D &);
  ObsVec1D & operator*= (const ObsVec1D &);
  ObsVec1D & operator/= (const ObsVec1D &);

  Eigen::VectorXd packEigen(const ObsVec1D &) const;
  size_t packEigenSize(const ObsVec1D &) const;

  size_t size() const {return data_.size();}
  const double & operator[](const std::size_t ii) const {return data_.at(ii);}
  double & operator[](const std::size_t ii) {return data_.at(ii);}

  void zero();
  /// set all values to ones (for tests)
  void ones();

  void axpy(const double &, const ObsVec1D &);
  void invert();
  void random();
  double dot_product_with(const ObsVec1D &) const;
  double rms() const;
  void mask(const ObsVec1D &);
  ObsVec1D & operator= (const ObsData1D<float> &);

  unsigned int nobs() const;
  const ObsTable & obsdb() const {return obsdb_;}

// I/O
  void save(const std::string &) const;
  void read(const std::string &);

  const double & missing() const {return missing_;}

 private:
  void print(std::ostream &) const;

  const ObsTable & obsdb_;
  std::vector<double> data_;
  const double missing_;
};
//-----------------------------------------------------------------------------
}  // namespace lorenz95
#endif  // LORENZ95_OBSVEC1D_H_
