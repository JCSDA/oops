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

#ifndef OOPS_ASSIMILATION_CONTROLVARIABLE_H_
#define OOPS_ASSIMILATION_CONTROLVARIABLE_H_

#include <cmath>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/Geometry.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/State.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace oops {

/// Control variable
/*!
 * The control variable acts as a container for the inputs of the variational
 * data assimilation cost functions in physical space.
 * That includes the states at the start the assimilation window or of each
 * sub-window but also additional variables such as model bias, VarBC
 * coefficients, or other control variables for algorithms that use them.
 * This is mostly a convenience class that is used to keep things together
 * and reduce the number of arguments to be passed around.
 */

template<typename MODEL, typename OBS> class ControlVariable;

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
class ControlVariable : public util::Printable,
                        public util::Serializable,
                        private util::ObjectCounter<ControlVariable<MODEL, OBS> > {
  typedef Geometry<MODEL>            Geometry_;
  typedef State4D<MODEL>             State4D_;
  typedef ModelAuxControl<MODEL>     ModelAux_;
  typedef ObsAuxControls<OBS>        ObsAux_;
  typedef ObsSpaces<OBS>             ObsSpaces_;

 public:
  static const std::string classname() {return "oops::ControlVariable";}

/// The arguments define the number of sub-windows and the resolution
  ControlVariable(std::shared_ptr<State4D_>,
                  std::shared_ptr<ModelAux_>, std::shared_ptr<ObsAux_>);
  explicit ControlVariable(const ControlVariable &, const bool copy = true);
  ~ControlVariable();

/// Get state control variable
  State<MODEL> & state(const size_t ii = 0) {return (*state_)[ii];}
  const State<MODEL> & state(const size_t ii = 0) const {return (*state_)[ii];}
  State4D<MODEL> & states() {return *state_;}
  const State4D<MODEL> & states() const {return *state_;}

/// Get augmented model control variable
  ModelAux_ & modVar() {return *modbias_;}
  const ModelAux_ & modVar() const {return *modbias_;}

/// Get augmented observation control variable
  ObsAux_ & obsVar() {return *obsbias_;}
  const ObsAux_ & obsVar() const {return *obsbias_;}

/// Serialize and deserialize ControlVariable
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

/// Assignment
  ControlVariable & operator= (const ControlVariable &);

 private:
  void print(std::ostream &) const override;

  std::shared_ptr<State4D_> state_;
  std::shared_ptr<ModelAux_> modbias_;  // not only for bias, better name?
  std::shared_ptr<ObsAux_> obsbias_;    // not only for bias, better name?
};

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
ControlVariable<MODEL, OBS>::ControlVariable(std::shared_ptr<State4D_> bg,
                                             std::shared_ptr<ModelAux_> maux,
                                             std::shared_ptr<ObsAux_> oaux)
  : state_(bg), modbias_(maux), obsbias_(oaux)
{
  Log::trace() << "ControlVariable::ControlVariable done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
ControlVariable<MODEL, OBS>::ControlVariable(const ControlVariable & other, const bool copy)
  : state_(new State4D_(*other.state_)), modbias_(), obsbias_()
{
  if (!copy) state_->zero();
  if (other.modbias_) modbias_.reset(new ModelAux_(*other.modbias_, copy));
  if (other.obsbias_) obsbias_.reset(new ObsAux_(*other.obsbias_, copy));
  Log::trace() << "ControlVariable copied" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
ControlVariable<MODEL, OBS>::~ControlVariable() {
  Log::trace() << "ControlVariable destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void ControlVariable<MODEL, OBS>::print(std::ostream & outs) const {
  outs << *state_ << std::endl;
  outs << *modbias_ << std::endl;
  outs << *obsbias_;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
size_t ControlVariable<MODEL, OBS>::serialSize() const {
  size_t ss = 0;
  for (size_t js = 0; js < state_->size(); ++js) {
    ss += (*state_)[js].serialSize();
  }
  ss += modbias_->serialSize();
  ss += obsbias_->serialSize();
  return ss;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void ControlVariable<MODEL, OBS>::serialize(std::vector<double> & vec) const {
  vec.reserve(vec.size() + this->serialSize());  // allocate memory to avoid reallocations
  for (size_t js = 0; js < state_->size(); ++js) {
    (*state_)[js].serialize(vec);
  }
  modbias_->serialize(vec);
  obsbias_->serialize(vec);
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void ControlVariable<MODEL, OBS>::deserialize(const std::vector<double> & vec, size_t & indx) {
  for (size_t js = 0; js < state_->size(); ++js) {
    (*state_)[js].deserialize(vec, indx);
  }
  modbias_->deserialize(vec, indx);
  obsbias_->deserialize(vec, indx);
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS> ControlVariable<MODEL, OBS> &
ControlVariable<MODEL, OBS>::operator=(const ControlVariable & rhs) {
    state_ = rhs.state_;
    modbias_ = rhs.modbias_;
    obsbias_ = rhs.obsbias_;
  return *this;
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_ASSIMILATION_CONTROLVARIABLE_H_
