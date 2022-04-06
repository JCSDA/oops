/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_CONTROLINCREMENT_H_
#define OOPS_ASSIMILATION_CONTROLINCREMENT_H_

#include <cmath>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/ObsAuxIncrements.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

namespace oops {

/// Control variable increment
/*!
 * The control variable acts as a container for the inputs of the variational
 * data assimilation cost functions in physical space.
 * That includes the states at the start the assimilation window or of each
 * sub-window but also additional variables such as model bias, VarBC
 * coefficients, or other control variables for algorithms that use them.
 * The control variable increment contains variations of the
 * control variable.
 */

template<typename MODEL, typename OBS> class CostJbTotal;

template<typename MODEL, typename OBS> class ControlIncrement;

// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
class ControlIncrement : public util::Printable,
                         public util::Serializable,
                         private util::ObjectCounter<ControlIncrement<MODEL, OBS> > {
  typedef CostJbTotal<MODEL, OBS>  JbTotal_;
  typedef Geometry<MODEL>          Geometry_;
  typedef Increment<MODEL>         Increment_;
  typedef ModelAuxIncrement<MODEL> ModelAuxIncr_;
  typedef ObsAuxIncrements<OBS>    ObsAuxIncrs_;

 public:
  static const std::string classname() {return "oops::ControlIncrement";}

/// Constructor, destructor
  explicit ControlIncrement(const JbTotal_ &);
  ControlIncrement(const ControlIncrement &, const bool copy = true);
  ControlIncrement(const ControlIncrement &, const eckit::Configuration &);
  ControlIncrement(const Geometry_ &, const ControlIncrement &);
  ~ControlIncrement();

/// Linear algebra operators
  void zero();
  ControlIncrement & operator=(const ControlIncrement &);
  ControlIncrement & operator+=(const ControlIncrement &);
  ControlIncrement & operator-=(const ControlIncrement &);
  ControlIncrement & operator*=(const double);
  void axpy(const double, const ControlIncrement &);
  double dot_product_with(const ControlIncrement &) const;

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;

/// Get geometry
  const Geometry_ & geometry() const {return increment_.geometry();}

/// Get state control variable
  Increment_ & state() {return increment_;}
  const Increment_ & state() const {return increment_;}

/// Get augmented model control variable
  ModelAuxIncr_ & modVar() {return modbias_;}
  const ModelAuxIncr_ & modVar() const {return modbias_;}

/// Get augmented observation control variable
  ObsAuxIncrs_ & obsVar() {return obsbias_;}
  const ObsAuxIncrs_ & obsVar() const {return obsbias_;}

/// Serialize and deserialize ControlIncrement
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

  void shift_forward();
  void shift_backward();

 private:
  void print(std::ostream &) const override;

  Increment_  increment_;
  ModelAuxIncr_ modbias_;   // not only for bias, better name?
  ObsAuxIncrs_  obsbias_;   // not only for bias, better name?
  const util::DateTime windowBegin_;
  const util::DateTime windowEnd_;
};

// =============================================================================

template<typename MODEL, typename OBS>
ControlIncrement<MODEL, OBS>::ControlIncrement(const JbTotal_ & jb)
  : increment_(*jb.jbState().newStateIncrement()),  // not good, extra copy
    modbias_(jb.resolution(), jb.jbModBias().config()),
    obsbias_(jb.jbObsBias().obspaces(), jb.jbObsBias().config()),
    windowBegin_(jb.windowBegin()), windowEnd_(jb.windowEnd())
{
  this->setObjectSize(this->serialSize()*sizeof(double));
  Log::trace() << "ControlIncrement:ControlIncrement created." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
ControlIncrement<MODEL, OBS>::ControlIncrement(const ControlIncrement & other, const bool copy)
  : increment_(other.increment_, copy), modbias_(other.modbias_, copy),
    obsbias_(other.obsbias_, copy), windowBegin_(other.windowBegin_), windowEnd_(other.windowEnd_)
{
  this->setObjectSize(this->serialSize()*sizeof(double));
  Log::trace() << "ControlIncrement:ControlIncrement copied." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
ControlIncrement<MODEL, OBS>::ControlIncrement(const ControlIncrement & other,
                                               const eckit::Configuration & tlConf)
  : increment_(other.increment_, tlConf), modbias_(other.modbias_, tlConf),
    obsbias_(other.obsbias_, tlConf), windowBegin_(other.windowBegin_), windowEnd_(other.windowEnd_)
{
  this->setObjectSize(this->serialSize()*sizeof(double));
  Log::trace() << "ControlIncrement:ControlIncrement copied." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
ControlIncrement<MODEL, OBS>::ControlIncrement(const Geometry_ & geom,
                                               const ControlIncrement & other)
  : increment_(geom, other.increment_), modbias_(other.modbias_, true),
    obsbias_(other.obsbias_, true), windowBegin_(other.windowBegin_), windowEnd_(other.windowEnd_)
{
  this->setObjectSize(this->serialSize()*sizeof(double));
  Log::trace() << "ControlIncrement:ControlIncrement copied." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
ControlIncrement<MODEL, OBS>::~ControlIncrement() {}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS> ControlIncrement<MODEL, OBS> &
ControlIncrement<MODEL, OBS>::operator=(const ControlIncrement & rhs) {
  increment_ = rhs.increment_;
  modbias_ = rhs.modbias_;
  obsbias_ = rhs.obsbias_;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS> ControlIncrement<MODEL, OBS> &
ControlIncrement<MODEL, OBS>::operator+=(const ControlIncrement & rhs) {
  increment_ += rhs.increment_;
  modbias_ += rhs.modbias_;
  obsbias_ += rhs.obsbias_;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS> ControlIncrement<MODEL, OBS> &
ControlIncrement<MODEL, OBS>::operator-=(const ControlIncrement & rhs) {
  increment_ -= rhs.increment_;
  modbias_ -= rhs.modbias_;
  obsbias_ -= rhs.obsbias_;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
ControlIncrement<MODEL, OBS> & ControlIncrement<MODEL, OBS>::operator*=(const double zz) {
  increment_ *= zz;
  modbias_ *= zz;
  obsbias_ *= zz;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
void ControlIncrement<MODEL, OBS>::zero() {
  increment_.zero();
  modbias_.zero();
  obsbias_.zero();
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
void ControlIncrement<MODEL, OBS>::axpy(const double zz, const ControlIncrement & rhs) {
  increment_.axpy(zz, rhs.increment_);
  modbias_.axpy(zz, rhs.modbias_);
  obsbias_.axpy(zz, rhs.obsbias_);
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
void ControlIncrement<MODEL, OBS>::read(const eckit::Configuration & config) {
  increment_.read(config);
  modbias_.read(config);
  obsbias_.read(config);
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
void ControlIncrement<MODEL, OBS>::write(const eckit::Configuration & config) const {
  increment_.write(config);
  modbias_.write(config);
  obsbias_.write(config);
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ControlIncrement<MODEL, OBS>::print(std::ostream & outs) const {
  outs << increment_;
  outs << modbias_;
  outs << obsbias_;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
double ControlIncrement<MODEL, OBS>::dot_product_with(const ControlIncrement & x2) const {
  double zz = 0.0;
  zz += dot_product(increment_, x2.increment_);
  zz += dot_product(modbias_, x2.modbias_);
  zz += dot_product(obsbias_, x2.obsbias_);
  return zz;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
size_t ControlIncrement<MODEL, OBS>::serialSize() const {
  size_t ss = 4;
  ss += increment_.serialSize();
  ss += modbias_.serialSize();
  ss += obsbias_.serialSize();
  return ss;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
void ControlIncrement<MODEL, OBS>::serialize(std::vector<double> & vec) const {
  vec.reserve(vec.size() + this->serialSize());  // allocate memory to avoid reallocations

  vec.push_back(-111.0);
  increment_.serialize(vec);

  vec.push_back(-222.0);
  modbias_.serialize(vec);

  vec.push_back(-333.0);
  obsbias_.serialize(vec);

  vec.push_back(-444.0);
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
void ControlIncrement<MODEL, OBS>::deserialize(const std::vector<double> & vec, size_t & indx) {
  ASSERT(vec.at(indx) == -111.0);
  ++indx;

  increment_.deserialize(vec, indx);

  ASSERT(vec.at(indx) == -222.0);
  ++indx;

  modbias_.deserialize(vec, indx);

  ASSERT(vec.at(indx) == -333.0);
  ++indx;

  obsbias_.deserialize(vec, indx);
  ASSERT(vec.at(indx) == -444.0);
  ++indx;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
void ControlIncrement<MODEL, OBS>::shift_forward() {
  increment_.shift_forward(windowBegin_);
// Probably needs some gathering of contributions for modbias_ and obsbias_
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
void ControlIncrement<MODEL, OBS>::shift_backward() {
  increment_.shift_backward(windowEnd_);
// Probably needs some gathering of contributions for modbias_ and obsbias_
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_CONTROLINCREMENT_H_
