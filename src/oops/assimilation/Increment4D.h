/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_INCREMENT4D_H_
#define OOPS_ASSIMILATION_INCREMENT4D_H_

#include <cmath>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/ptr_container/ptr_map.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/CostJbState.h"
#include "oops/assimilation/State4D.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

/// State increment
/*!
 * The state increment contains the increment to the 3D or 4D state part of
 * the VDA control variable.
 */

// -----------------------------------------------------------------------------
template<typename MODEL> class Increment4D : public util::Printable {
  typedef Geometry<MODEL>    Geometry_;
  typedef Increment<MODEL>   Increment_;
  typedef CostJbState<MODEL> JbState_;
  typedef State4D<MODEL>     State4D_;

 public:
  static const std::string classname() {return "Increment4D";}

/// Constructor, destructor
  explicit Increment4D(const JbState_ &);
  Increment4D(const Increment4D &, const bool copy = true);
  Increment4D(const Increment4D &, const eckit::Configuration &);
  Increment4D(const Geometry_ &, const Increment4D &);
  ~Increment4D();

/// Linear algebra operators
  void diff(const State4D_ &, const State4D_ &);
  void zero();
  Increment4D & operator=(const Increment4D &);
  Increment4D & operator+=(const Increment4D &);
  Increment4D & operator-=(const Increment4D &);
  Increment4D & operator*=(const double);
  void axpy(const double, const Increment4D &);
  double dot_product_with(const Increment4D &) const;

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;

/// Get geometry
  Geometry_ geometry() const {return this->get(first_).geometry();}

/// Get model space control variable
  Increment_ & operator[](const int ii) {return this->get(ii);}
  const Increment_ & operator[](const int ii) const {return this->get(ii);}
  int first() const {return first_;}
  int last() const {return last_;}

/// To be removed
  void shift_forward();
  void shift_backward();

/// Serialize-Deserialize an Increment4D
  void serialize(std::vector<double> &) const;
  void deserialize(const std::vector<double> &);

 private:
  Increment_ & get(const int);
  const Increment_ & get(const int) const;
  void print(std::ostream &) const;

  typedef typename boost::ptr_map<int, Increment_>::iterator iter_;
  typedef typename boost::ptr_map<int, Increment_>::const_iterator icst_;
  boost::ptr_map<int, Increment_> incr4d_;
  int first_;
  int last_;
};

// =============================================================================

template<typename MODEL>
Increment4D<MODEL>::Increment4D(const JbState_ & jb)
  : incr4d_(), first_(0), last_(jb.nstates() - 1)
{
  for (int jsub = 0; jsub <= last_; ++jsub) {
    Increment_ * incr = jb.newStateIncrement(jsub);
    incr4d_.insert(jsub, incr);
  }
  Log::trace() << "Increment4D:Increment4D created." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Increment4D<MODEL>::Increment4D(const Increment4D & other, const bool copy)
  : incr4d_(), first_(other.first_), last_(other.last_)
{
  for (icst_ jsub = other.incr4d_.begin(); jsub != other.incr4d_.end(); ++jsub) {
    int isub = jsub->first;
    Increment_ * tmp = new Increment_(*jsub->second, copy);
    incr4d_.insert(isub, tmp);
  }
  Log::trace() << "Increment4D:Increment4D copied." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Increment4D<MODEL>::Increment4D(const Increment4D & other, const eckit::Configuration & tlConf)
  : incr4d_(), first_(other.first_), last_(other.last_)
{
  for (icst_ jsub = other.incr4d_.begin(); jsub != other.incr4d_.end(); ++jsub) {
    int isub = jsub->first;
    Increment_ * tmp = new Increment_(*jsub->second);
    incr4d_.insert(isub, tmp);
  }
  Log::trace() << "Increment4D:Increment4D copied." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Increment4D<MODEL>::Increment4D(const Geometry_ & geom, const Increment4D & other)
  : incr4d_(), first_(other.first_), last_(other.last_)
{
  for (icst_ jsub = other.incr4d_.begin(); jsub != other.incr4d_.end(); ++jsub) {
    int isub = jsub->first;
    Increment_ * tmp = new Increment_(geom, *jsub->second);
    incr4d_.insert(isub, tmp);
  }
  Log::trace() << "Increment4D:Increment4D copied." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Increment<MODEL> & Increment4D<MODEL>::get(const int ii) {
  iter_ it = incr4d_.find(ii);
  ASSERT(it != incr4d_.end());
  return *it->second;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
const Increment<MODEL> & Increment4D<MODEL>::get(const int ii) const {
  icst_ it = incr4d_.find(ii);
  ASSERT(it != incr4d_.end());
  return *it->second;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Increment4D<MODEL>::~Increment4D() {}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::zero() {
  for (iter_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    jsub->second->zero();
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::diff(const State4D_ & cv1, const State4D_ & cv2) {
  for (iter_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    jsub->second->diff(cv1[jsub->first], cv2[jsub->first]);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Increment4D<MODEL> & Increment4D<MODEL>::operator=(const Increment4D & rhs) {
  incr4d_ = rhs.incr4d_;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Increment4D<MODEL> & Increment4D<MODEL>::operator+=(const Increment4D & rhs) {
  for (int jsub = rhs.first(); jsub <= rhs.last(); ++jsub) {
    this->get(jsub) += rhs[jsub];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Increment4D<MODEL> & Increment4D<MODEL>::operator-=(const Increment4D & rhs) {
  for (int jsub = rhs.first(); jsub <= rhs.last(); ++jsub) {
    this->get(jsub) -= rhs[jsub];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Increment4D<MODEL> & Increment4D<MODEL>::operator*=(const double zz) {
  for (iter_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    *jsub->second *= zz;
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::axpy(const double zz, const Increment4D & rhs) {
  for (iter_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    jsub->second->axpy(zz, rhs[jsub->first]);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::read(const eckit::Configuration & config) {
  for (iter_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    std::stringstream ss;
    ss << jsub->first+1;
    std::string query = "increment[@indx='" + ss.str() + "']";
    eckit::LocalConfiguration fileConfig(config, query);
    jsub->second->read(fileConfig);
    Log::info() << "Increment4D:read increment" << *jsub->second << std::endl;
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::write(const eckit::Configuration & config) const {
  for (icst_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    jsub->second->write(config);
  }
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Increment4D<MODEL>::print(std::ostream & outs) const {
  for (icst_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    outs << *jsub->second << std::endl;
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
double Increment4D<MODEL>::dot_product_with(const Increment4D & x2) const {
  double zz = 0.0;
  for (icst_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    zz += dot_product(*jsub->second, x2[jsub->first]);
  }
  return zz;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::shift_forward() {
  typedef typename boost::ptr_map<int, Increment_>::reverse_iterator rit;
  for (rit jsub = incr4d_.rbegin(); jsub != incr4d_.rend(); ++jsub) {
    const int isub = jsub->first;
    if (isub > first_) this->get(isub) = this->get(isub-1);
  }
  incr4d_.erase(first_);
  Log::info() << "Increment4D::shift_forward erased " << first_ << std::endl;
  first_ += 1;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::shift_backward() {
  for (iter_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    const int isub = jsub->first;
    if (isub < last_) this->get(isub) = this->get(isub+1);
  }
  incr4d_.erase(last_);
  Log::info() << "Increment4D::shift_backward erased " << last_ << std::endl;
  last_ -= 1;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::serialize(std::vector<double> & vect) const {
  for (iter_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    jsub->second->serialize(vect);
  }
  Log::info() << "Increment4D::serialize done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::deserialize(const std::vector<double> & vect) {
  int size_vec = 0;
  int sum = 0;
  int ii = 0;
  int bgn = 0;

  for (iter_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    sum += size_vec;
    bgn = sum + 1 + ii;  // bgn = 1  for the first increment
    size_vec = std::lround(vect[bgn - 1]);
    std::vector<double> incr_ii(vect.begin() + bgn, vect.begin() + bgn + size_vec);
    jsub->second->deserialize(incr_ii);
    ++ii;
  }
  Log::info() << "Increment4D::deserialize done" << std::endl;
}

// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_ASSIMILATION_INCREMENT4D_H_
