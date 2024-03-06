/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_OBSVECQG_H_
#define QG_MODEL_OBSVECQG_H_

#include <Eigen/Dense>
#include <ostream>
#include <string>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "oops/qg/QgFortran.h"

namespace qg {
  class ObsSpaceQG;
  template <typename DATATYPE> class ObsDataQG;

// -----------------------------------------------------------------------------
/// ObsVecQG class to handle vectors in observation space for QG model.

class ObsVecQG : public util::Printable,
                 private util::ObjectCounter<ObsVecQG> {
 public:
  static const std::string classname() {return "qg::ObsVecQG";}

  explicit ObsVecQG(const ObsSpaceQG &,
           const std::string & name = std::string());
  ObsVecQG(const ObsVecQG &);
  ~ObsVecQG();

  ObsVecQG & operator = (const ObsVecQG &);
  ObsVecQG & operator*= (const double &);
  ObsVecQG & operator+= (const ObsVecQG &);
  ObsVecQG & operator-= (const ObsVecQG &);
  ObsVecQG & operator*= (const ObsVecQG &);
  ObsVecQG & operator/= (const ObsVecQG &);

  Eigen::VectorXd packEigen(const ObsVecQG &) const;
  size_t packEigenSize(const ObsVecQG &) const;
  size_t size() const;

  /// set all values to zero
  void zero();
  /// set \p i-th value to missing value
  void setToMissing(int i);
  /// set all values to one
  void ones();
  void axpy(const double &, const ObsVecQG &);
  void invert();
  void random();
  double dot_product_with(const ObsVecQG &) const;
  double rms() const;
  void mask(const ObsVecQG &);
  ObsVecQG & operator=(const ObsDataQG<float> &);

  unsigned int nobs() const;

  const int & toFortran() const {return keyOvec_;}

// I/O
  void save(const std::string &) const;
  void read(const std::string &);

 private:
  void print(std::ostream &) const;

  const ObsSpaceQG & obsdb_;
  F90ovec keyOvec_;
};
// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_OBSVECQG_H_
