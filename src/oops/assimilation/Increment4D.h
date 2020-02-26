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

#include "atlas/field/FieldSet.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/assimilation/CostJbState.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/Variables.h"
#if !ATLASIFIED
#include "oops/generic/UnstructuredGrid.h"
#endif
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
  explicit Increment4D(const Increment_ &);
  Increment4D(const Geometry_ &, const Variables &, const std::vector<util::DateTime> &);
  Increment4D(const Increment4D &, const bool copy = true);
  Increment4D(const Geometry_ &, const Increment4D &);
  ~Increment4D();

/// Interfacing
  Increment_ & incr4d() {return *incr4d_;}
  const Increment_ & incr4d() const {return *incr4d_;}

/// Linear algebra operators
  void diff(const State4D_ &, const State4D_ &);
  void zero();
  void dirac(std::vector<eckit::LocalConfiguration>);
  Increment4D & operator=(const Increment4D &);
  Increment4D & operator+=(const Increment4D &);
  Increment4D & operator-=(const Increment4D &);
  Increment4D & operator*=(const double);
  void axpy(const double, const Increment4D &, const bool check = true);
  double dot_product_with(const Increment4D &) const;
  void schur_product_with(const Increment4D &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;

/// Get geometry
  Geometry_ geometry() const {return this->get(first_).geometry();}

#if ATLASIFIED
/// ATLAS FieldSet
  void setAtlas(atlas::FieldSet *) const;
  void toAtlas(atlas::FieldSet *) const;
  void fromAtlas(atlas::FieldSet *);
#else
/// Unstructured grid
  void ug_coord(UnstructuredGrid &) const;
  void field_to_ug(UnstructuredGrid &) const;
  void field_from_ug(const UnstructuredGrid &);
#endif

/// Get model space control variable
  Increment_ & operator[](const int ii) {return this->get(ii);}
  const Increment_ & operator[](const int ii) const {return this->get(ii);}
  int first() const {return first_;}
  int last() const {return last_;}
  size_t size() const {return last_-first_+1;}

/// To be removed
  void shift_forward();
  void shift_backward();

/// Serialize and deserialize
  size_t serialSize() const;
  void serialize(std::vector<double> &) const;
  void deserialize(const std::vector<double> &, size_t &);

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
template <typename MODEL>
State4D<MODEL> & operator+=(State4D<MODEL> & xx, const Increment4D<MODEL> & dx) {
  Log::trace() << "operator+=(State4D, Increment4D) starting" << std::endl;
  for (size_t ii = 0; ii < xx.size(); ++ii) {
    xx[ii] += dx[ii];
  }
  Log::trace() << "operator+=(State4D, Increment4D) done" << std::endl;
  return xx;
}
// ----------------------------------------------------------------------------
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
Increment4D<MODEL>::Increment4D(const Increment_ & dx)
  : incr4d_(), first_(0), last_(0)
{
  Increment_ * incr = new Increment_(dx);
  incr4d_.insert(0, incr);
  Log::trace() << "Increment4D:Increment4D created." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Increment4D<MODEL>::Increment4D(const Geometry_ & resol,
                                const Variables & vars,
                                const std::vector<util::DateTime> & timeslots)
  : incr4d_(), first_(0), last_(timeslots.size() - 1)
{
  for (int jsub = 0; jsub <= last_; ++jsub) {
    Increment_ * incr = new Increment_(resol, vars, timeslots[jsub]);
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
void Increment4D<MODEL>::dirac(std::vector<eckit::LocalConfiguration> confs) {
  this->zero();
  for (const auto & conf : confs) {
    const util::DateTime date(conf.getString("date"));
    for (iter_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
      if (date == jsub->second->validTime()) {
        jsub->second->dirac(conf);
      }
    }
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
void Increment4D<MODEL>::axpy(const double zz, const Increment4D & rhs, const bool check) {
  for (iter_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    jsub->second->axpy(zz, rhs[jsub->first], check);
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
#if ATLASIFIED
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::setAtlas(atlas::FieldSet * afieldset) const {
  for (icst_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    jsub->second->setAtlas(afieldset);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::toAtlas(atlas::FieldSet * afieldset) const {
  for (icst_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    jsub->second->toAtlas(afieldset);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::fromAtlas(atlas::FieldSet * afieldset) {
  for (iter_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    jsub->second->fromAtlas(afieldset);
  }
}
// -----------------------------------------------------------------------------
#else
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::ug_coord(UnstructuredGrid & ug) const {
  incr4d_.begin()->second->increment().ug_coord(ug);
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::field_to_ug(UnstructuredGrid & ug) const {
  for (icst_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    jsub->second->increment().field_to_ug(ug, jsub->first);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::field_from_ug(const UnstructuredGrid & ug) {
  for (iter_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    jsub->second->increment().field_from_ug(ug, jsub->first);
  }
}
// -----------------------------------------------------------------------------
#endif
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
void Increment4D<MODEL>::schur_product_with(const Increment4D & x2) {
  for (iter_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    jsub->second->schur_product_with(x2[jsub->first]);
  }
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
size_t Increment4D<MODEL>::serialSize() const {
  size_t ss = 1;
  for (icst_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    ss += jsub->second->serialSize();
    ++ss;
  }
  return ss;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::serialize(std::vector<double> & vect) const {
  vect.push_back(-98765.4321);
  for (icst_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    jsub->second->serialize(vect);
    vect.push_back(-98765.4321);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment4D<MODEL>::deserialize(const std::vector<double> & vect, size_t & current) {
  ASSERT(vect.at(current) == -98765.4321);
  ++current;
  for (iter_ jsub = incr4d_.begin(); jsub != incr4d_.end(); ++jsub) {
    jsub->second->deserialize(vect, current);
    ASSERT(vect.at(current) == -98765.4321);
    ++current;
  }
}

// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_ASSIMILATION_INCREMENT4D_H_
