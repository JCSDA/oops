/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_INCREMENT4D_H_
#define OOPS_BASE_INCREMENT4D_H_

#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

/// 4D model state Increment (vector of 3D Increments)
template<typename MODEL> class Increment4D : public util::Printable {
  typedef Geometry<MODEL>    Geometry_;
  typedef Increment<MODEL>   Increment_;
  typedef State4D<MODEL>     State4D_;

 public:
  static const std::string classname() {return "Increment4D";}

  /// Constructor for specified times
  Increment4D(const Geometry_ &, const Variables &, const std::vector<util::DateTime> &);

  /// Linear algebra operators
  void diff(const State4D_ &, const State4D_ &);
  void zero();
  void random();
  void ones();
  double dot_product_with(const Increment4D &) const;
  void schur_product_with(const Increment4D &);

  /// Get geometry
  Geometry_ geometry() const {return incr4d_[0].geometry();}

  /// Get 3D increments
  Increment_ & operator[](const int ii) {return incr4d_[ii];}
  const Increment_ & operator[](const int ii) const {return incr4d_[ii];}
  size_t size() const {return incr4d_.size();}

 private:
  void print(std::ostream &) const override;

  std::vector<Increment_> incr4d_;
};

// =============================================================================
/// "Increment" 4D State \p xx with 4D Increment \p dx
template <typename MODEL>
State4D<MODEL> & operator+=(State4D<MODEL> & xx, const Increment4D<MODEL> & dx) {
  Log::trace() << "operator+=(State4D, Increment4D) starting" << std::endl;
  for (size_t ii = 0; ii < xx.size(); ++ii) {
    xx[ii] += dx[ii];
  }
  Log::trace() << "operator+=(State4D, Increment4D) done" << std::endl;
  return xx;
}

// -----------------------------------------------------------------------------
template<typename MODEL>
Increment4D<MODEL>::Increment4D(const Geometry_ & resol,
                                const Variables & vars,
                                const std::vector<util::DateTime> & timeslots)
  : incr4d_()
{
  for (size_t jtime = 0; jtime < timeslots.size(); ++jtime) {
    incr4d_.emplace_back(resol, vars, timeslots[jtime]);
  }
  Log::trace() << "Increment4D:Increment4D created." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::zero() {
  for (auto & incr : incr4d_) {
    incr.zero();
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::random() {
  for (auto & incr : incr4d_) {
    incr.random();
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::ones() {
  for (auto & incr : incr4d_) {
    incr.ones();
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::diff(const State4D_ & cv1, const State4D_ & cv2) {
  for (size_t jtime = 0; jtime < incr4d_.size(); ++jtime) {
    incr4d_[jtime].diff(cv1[jtime], cv2[jtime]);
  }
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Increment4D<MODEL>::print(std::ostream & outs) const {
  for (const auto & incr : incr4d_) {
    outs << incr << std::endl;
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
double Increment4D<MODEL>::dot_product_with(const Increment4D & x2) const {
  double zz = 0.0;
  for (size_t jtime = 0; jtime < incr4d_.size(); ++jtime) {
    zz += dot_product(incr4d_[jtime], x2[jtime]);
  }
  return zz;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::schur_product_with(const Increment4D & x2) {
  for (size_t jtime = 0; jtime < incr4d_.size(); ++jtime) {
    incr4d_[jtime].schur_product_with(x2[jtime]);
  }
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_INCREMENT4D_H_
