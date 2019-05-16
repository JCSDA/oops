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
#include "oops/assimilation/Increment4D.h"
#include "oops/base/ObsAuxIncrements.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

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

template<typename MODEL> class CostJbTotal;

template<typename MODEL> class ControlIncrement;

// -----------------------------------------------------------------------------
template<typename MODEL>
class ControlIncrement : public util::Printable,
                         private util::ObjectCounter<ControlIncrement<MODEL> > {
  typedef CostJbTotal<MODEL>       JbTotal_;
  typedef Geometry<MODEL>          Geometry_;
  typedef Increment4D<MODEL>       Increment4D_;
  typedef ModelAuxIncrement<MODEL> ModelAuxIncr_;
  typedef ObsAuxIncrements<MODEL>  ObsAuxIncrs_;

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
  Geometry_ geometry() const {return incrm4d_.geometry();}

/// Get state control variable
  Increment4D_ & state() {return incrm4d_;}
  const Increment4D_ & state() const {return incrm4d_;}

/// Get augmented model control variable
  ModelAuxIncr_ & modVar() {return modbias_;}
  const ModelAuxIncr_ & modVar() const {return modbias_;}

/// Get augmented observation control variable
  ObsAuxIncrs_ & obsVar() {return obsbias_;}
  const ObsAuxIncrs_ & obsVar() const {return obsbias_;}

/// Serialize and deserialize ControlIncrement
  std::vector<double> serialize() const;
  void deserialize(const std::vector<double> &);

 private:
  void print(std::ostream &) const;
  Increment4D_  incrm4d_;
  ModelAuxIncr_ modbias_;   // not only for bias, better name?
  ObsAuxIncrs_  obsbias_;   // not only for bias, better name?
};

// =============================================================================

template<typename MODEL>
ControlIncrement<MODEL>::ControlIncrement(const JbTotal_ & jb)
  : incrm4d_(jb.jbState()), modbias_(jb.resolution(), jb.jbModBias().config()),
    obsbias_(jb.jbObsBias().config())
{
  Log::trace() << "ControlIncrement:ControlIncrement created." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ControlIncrement<MODEL>::ControlIncrement(const ControlIncrement & other, const bool copy)
  : incrm4d_(other.incrm4d_, copy), modbias_(other.modbias_, copy),
    obsbias_(other.obsbias_, copy)
{
  Log::trace() << "ControlIncrement:ControlIncrement copied." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ControlIncrement<MODEL>::ControlIncrement(const ControlIncrement & other,
                                          const eckit::Configuration & tlConf)
  : incrm4d_(other.incrm4d_, tlConf), modbias_(other.modbias_, tlConf),
    obsbias_(other.obsbias_, tlConf)
{
  Log::trace() << "ControlIncrement:ControlIncrement copied." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ControlIncrement<MODEL>::ControlIncrement(const Geometry_ & geom,
                                          const ControlIncrement & other)
  : incrm4d_(geom, other.incrm4d_), modbias_(other.modbias_, true),
    obsbias_(other.obsbias_, true)
{
  Log::trace() << "ControlIncrement:ControlIncrement copied." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ControlIncrement<MODEL>::~ControlIncrement() {}
// -----------------------------------------------------------------------------
template<typename MODEL> ControlIncrement<MODEL> &
ControlIncrement<MODEL>::operator=(const ControlIncrement & rhs) {
  incrm4d_ = rhs.incrm4d_;
  modbias_ = rhs.modbias_;
  obsbias_ = rhs.obsbias_;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL> ControlIncrement<MODEL> &
ControlIncrement<MODEL>::operator+=(const ControlIncrement & rhs) {
  incrm4d_ += rhs.incrm4d_;
  modbias_ += rhs.modbias_;
  obsbias_ += rhs.obsbias_;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL> ControlIncrement<MODEL> &
ControlIncrement<MODEL>::operator-=(const ControlIncrement & rhs) {
  incrm4d_ -= rhs.incrm4d_;
  modbias_ -= rhs.modbias_;
  obsbias_ -= rhs.obsbias_;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ControlIncrement<MODEL> & ControlIncrement<MODEL>::operator*=(const double zz) {
  incrm4d_ *= zz;
  modbias_ *= zz;
  obsbias_ *= zz;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ControlIncrement<MODEL>::zero() {
  incrm4d_.zero();
  modbias_.zero();
  obsbias_.zero();
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ControlIncrement<MODEL>::axpy(const double zz, const ControlIncrement & rhs) {
  incrm4d_.axpy(zz, rhs.incrm4d_);
  modbias_.axpy(zz, rhs.modbias_);
  obsbias_.axpy(zz, rhs.obsbias_);
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ControlIncrement<MODEL>::read(const eckit::Configuration & config) {
  incrm4d_.read(config);
  modbias_.read(config);
  obsbias_.read(config);
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ControlIncrement<MODEL>::write(const eckit::Configuration & config) const {
  incrm4d_.write(config);
  modbias_.write(config);
  obsbias_.write(config);
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ControlIncrement<MODEL>::print(std::ostream & outs) const {
  outs << incrm4d_;
  outs << modbias_;
  outs << obsbias_;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
double ControlIncrement<MODEL>::dot_product_with(const ControlIncrement & x2) const {
  double zz = 0.0;
  zz += dot_product(incrm4d_, x2.incrm4d_);
  zz += dot_product(modbias_, x2.modbias_);
  zz += dot_product(obsbias_, x2.obsbias_);
  return zz;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
std::vector<double> ControlIncrement<MODEL>::serialize() const {
  std::vector<double> vec;
  incrm4d_.serialize(vec);
  vec.insert(vec.begin(), static_cast<double>(vec.size()));  // includes info+incr+time and date
  vec.push_back(500.00);
  modbias_.serialize(vec);
  vec.push_back(500.00);
  obsbias_.serialize(vec);
  return vec;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ControlIncrement<MODEL>::deserialize(const std::vector<double> & vec) {
  unsigned int s_incrm4d = std::lround(vec[0]);
  unsigned int s_modbias = std::lround(vec[s_incrm4d + 2]);
  unsigned int s_obsbias = std::lround(vec[s_incrm4d + s_modbias + 4]);

  ASSERT(vec[s_incrm4d + 1] == 500);
  ASSERT(vec[s_incrm4d + s_modbias + 3] == 500);

  std::vector<double> vec_incrm4d(vec.begin() + 1, vec.begin() + 1 + s_incrm4d);
  std::vector<double> vec_modbias(vec.begin() + s_incrm4d + 3,
                                  vec.begin() + s_incrm4d + 3 + s_modbias);
  std::vector<double> vec_obsbias(vec.begin() + s_incrm4d + s_modbias + 5,
                                  vec.begin() + s_incrm4d + s_modbias + 5 + s_obsbias);

  incrm4d_.deserialize(vec_incrm4d);
  modbias_.deserialize(vec_modbias);
  obsbias_.deserialize(vec_obsbias);
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_CONTROLINCREMENT_H_
