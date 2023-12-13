/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Crown Copyright 2023, the Met Office.
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
#include "oops/assimilation/ControlVariable.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/Increment4D.h"
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
  typedef CostJbTotal<MODEL, OBS>      JbTotal_;
  typedef ControlVariable<MODEL, OBS>  CtrlVar_;
  typedef Geometry<MODEL>              Geometry_;
  typedef Increment<MODEL>             Increment_;
  typedef ModelAuxIncrement<MODEL>     ModelAuxIncr_;
  typedef ObsAuxIncrements<OBS>        ObsAuxIncrs_;

 public:
  static const std::string classname() {return "oops::ControlIncrement";}

/// Constructor, destructor
  explicit ControlIncrement(const JbTotal_ &);
  ControlIncrement(const ControlIncrement &, const bool copy = true);
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
  void schur_product_with(const ControlIncrement & other);

  /// Set this ControlIncrement to be difference between \p cvar1 and \p cvar2
  void diff(const CtrlVar_ & cvar1, const CtrlVar_ & cvar2);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;

/// Get geometry
  const Geometry_ & geometry() const {return increment_.geometry();}

/// Get state control variable
  Increment_ & state(const size_t ii = 0) {return increment_[ii];}
  const Increment_ & state(const size_t ii = 0) const {return increment_[ii];}
  Increment4D<MODEL> & states() {return increment_;}
  const Increment4D<MODEL> & states() const {return increment_;}

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

  Increment4D<MODEL>  increment_;
  ModelAuxIncr_ modbias_;   // not only for bias, better name?
  ObsAuxIncrs_  obsbias_;   // not only for bias, better name?
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
ControlVariable<MODEL, OBS> & operator+=(ControlVariable<MODEL, OBS> & xx,
                                    const ControlIncrement<MODEL, OBS> & dx) {
  Log::trace() << "operator+=(ControlVariable, ControlIncrement) starting" << std::endl;
  util::Timer timer("oops::ControlIncrement", "operator+=ControlVariable");
  xx.state() += dx.state();
  xx.modVar() += dx.modVar();
  xx.obsVar() += dx.obsVar();
  Log::trace() << "operator+=(ControlVariable, ControlIncrement) done" << std::endl;
  return xx;
}

// =============================================================================

template<typename MODEL, typename OBS>
ControlIncrement<MODEL, OBS>::ControlIncrement(const JbTotal_ & jb)
  : increment_(jb.jbState().geometry(), jb.jbState().variables(), jb.jbState().times(),
                                                                  jb.jbState().comm()),
    modbias_(jb.jbModBias().geometry(), jb.jbModBias().config()),
    obsbias_(jb.jbObsBias().obspaces(), jb.jbObsBias().config())
{
  this->setObjectSize(this->serialSize()*sizeof(double));
  Log::trace() << "ControlIncrement:ControlIncrement created." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
ControlIncrement<MODEL, OBS>::ControlIncrement(const ControlIncrement & other, const bool copy)
  : increment_(other.increment_, copy),
    modbias_(other.modbias_, copy), obsbias_(other.obsbias_, copy)
{
  this->setObjectSize(this->serialSize()*sizeof(double));
  Log::trace() << "ControlIncrement:ControlIncrement copied." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
ControlIncrement<MODEL, OBS>::ControlIncrement(const Geometry_ & geom,
                                               const ControlIncrement & other)
  : increment_(geom, other.increment_),
    modbias_(other.modbias_, true), obsbias_(other.obsbias_, true)
{
  this->setObjectSize(this->serialSize()*sizeof(double));
  Log::trace() << "ControlIncrement:ControlIncrement copied." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
ControlIncrement<MODEL, OBS>::~ControlIncrement() {}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
void ControlIncrement<MODEL, OBS>::diff(const CtrlVar_ & cvar1, const CtrlVar_ & cvar2) {
  increment_.diff(cvar1.states(), cvar2.states());
  modbias_.diff(cvar1.modVar(), cvar2.modVar());
  obsbias_.diff(cvar1.obsVar(), cvar2.obsVar());
}
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
  ASSERT(increment_.is_3d());  // for now
  increment_[0].read(config);
  modbias_.read(config);
  obsbias_.read(config);
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
void ControlIncrement<MODEL, OBS>::write(const eckit::Configuration & config) const {
  increment_.write(config.getSubConfiguration("state component"));
  modbias_.write(config.getSubConfiguration("model aux component"));
  obsbias_.write(config.getSubConfiguration("obs aux component"));
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
void ControlIncrement<MODEL, OBS>::schur_product_with(const ControlIncrement & other) {
  increment_.schur_product_with(other.increment_);
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
size_t ControlIncrement<MODEL, OBS>::serialSize() const {
  size_t ss = 0;
  for (size_t js = 0; js < increment_.size(); ++js) {
    ss += increment_[js].serialSize();
  }
  ss += modbias_.serialSize();
  ss += obsbias_.serialSize();
  return ss;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
void ControlIncrement<MODEL, OBS>::serialize(std::vector<double> & vec) const {
  vec.reserve(vec.size() + this->serialSize());  // allocate memory to avoid reallocations
  for (size_t js = 0; js < increment_.size(); ++js) {
    increment_[js].serialize(vec);
  }
  modbias_.serialize(vec);
  obsbias_.serialize(vec);
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
void ControlIncrement<MODEL, OBS>::deserialize(const std::vector<double> & vec, size_t & indx) {
  for (size_t js = 0; js < increment_.size(); ++js) {
    increment_[js].deserialize(vec, indx);
  }
  modbias_.deserialize(vec, indx);
  obsbias_.deserialize(vec, indx);
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
void ControlIncrement<MODEL, OBS>::shift_forward() {
  increment_.shift_forward();
// Probably needs some gathering of contributions for modbias_ and obsbias_
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
void ControlIncrement<MODEL, OBS>::shift_backward() {
  increment_.shift_backward();
// Probably needs some gathering of contributions for modbias_ and obsbias_
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_CONTROLINCREMENT_H_
