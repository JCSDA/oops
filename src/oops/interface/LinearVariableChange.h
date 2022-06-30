/*
 * (C) Copyright 2018-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_LINEARVARIABLECHANGE_H_
#define OOPS_INTERFACE_LINEARVARIABLECHANGE_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/LinearVariableChangeParametersBase.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/GenericParameters.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParametersOrConfiguration.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// \brief MODEL-agnostic part of the linear variable change
///
/// Note: each implementation should typedef `LinVariableChangeParams_` to the name of a
/// subclass of oops::LinearVariableChangeParametersBase holding its configuration settings and
/// provide a constructor with the following signature:
///
///     LinearVariableChange(const Geometry_ &, const LinVariableChangeParams_ &);
template <typename MODEL>
class LinearVariableChange {
  typedef typename MODEL::LinearVariableChange LinearVariableChange_;
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
  /// A subclass of oops::LinearVariableChangeParametersBase holding the configuration settings of
  /// the variable change.
  /// Defined as LinearVariableChange_::Parameters_ if LinearVariableChange_ defines a Parameters_
  /// type; otherwise as GenericLinearVariableChangeParameters.
  typedef TParameters_IfAvailableElseFallbackType_t<LinearVariableChange_,
                                    GenericLinearVariableChangeParameters> Parameters_;
  static const std::string classname() {return "oops::LinearVariableChange";}

  LinearVariableChange(const Geometry_ &, const Parameters_ &);
  LinearVariableChange(const Geometry_ &, const eckit::Configuration &);

  virtual ~LinearVariableChange();

  void setTrajectory(const State_ &, const State_ &);
  void setTrajectory(const State_ &);
  void multiply(Increment_ &, const Variables &) const;
  void multiplyInverse(Increment_ &, const Variables &) const;
  void multiplyAD(Increment_ &, const Variables &) const;
  void multiplyInverseAD(Increment_ &, const Variables &) const;

 private:
  void print(std::ostream &) const;

  std::unique_ptr<LinearVariableChange_> chvar_;
};

// =============================================================================

template<typename MODEL>
LinearVariableChange<MODEL>::LinearVariableChange(const Geometry_ & resol,
    const eckit::Configuration & conf)
  :  LinearVariableChange<MODEL>::LinearVariableChange(resol,
        validateAndDeserialize<Parameters_>(conf))
{}

template<typename MODEL>
LinearVariableChange<MODEL>::LinearVariableChange(const Geometry_ & resol,
    const Parameters_ & parameters) : chvar_() {
  Log::trace() << "LinearVariableChange<MODEL>::LinearVariableChange starting" << std::endl;
  util::Timer timer(classname(), "LinearVariableChange");
  chvar_.reset(new LinearVariableChange_(resol.geometry(), parameters));
  Log::trace() << "LinearVariableChange<MODEL>::LinearVariableChange done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
LinearVariableChange<MODEL>::~LinearVariableChange() {
  Log::trace() << "LinearVariableChange<MODEL>::~LinearVariableChange starting" << std::endl;
  util::Timer timer(classname(), "~LinearVariableChange");
  chvar_.reset();
  Log::trace() << "LinearVariableChange<MODEL>::~LinearVariableChange done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearVariableChange<MODEL>::multiply(Increment_ & dx, const Variables & vars) const {
  Log::trace() << "LinearVariableChange<MODEL>::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  chvar_->multiply(dx.increment(), vars);
  Log::trace() << "LinearVariableChange<MODEL>::multiply done" << std::endl;
}
// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearVariableChange<MODEL>::multiplyInverse(Increment_ & dx,
                                                const Variables & vars) const {
  Log::trace() << "LinearVariableChange<MODEL>::multiplyInverse starting" << std::endl;
  util::Timer timer(classname(), "multiplyInverse");
  chvar_->multiplyInverse(dx.increment(), vars);
  Log::trace() << "LinearVariableChange<MODEL>::multiplyInverse done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearVariableChange<MODEL>::multiplyAD(Increment_ & dx,
                                                const Variables & vars) const {
  Log::trace() << "LinearVariableChange<MODEL>::multiplyAD starting" << std::endl;
  util::Timer timer(classname(), "multiplyAD");
  chvar_->multiplyAD(dx.increment(), vars);
  Log::trace() << "LinearVariableChange<MODEL>::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearVariableChange<MODEL>::multiplyInverseAD(Increment_ & dx,
                                                const Variables & vars) const {
  Log::trace() << "LinearVariableChange<MODEL>::multiplyInverseAD starting" << std::endl;
  util::Timer timer(classname(), "multiplyInverseAD");
  chvar_->multiplyInverseAD(dx.increment(), vars);
  Log::trace() << "LinearVariableChange<MODEL>::multiplyInverseAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearVariableChange<MODEL>::setTrajectory(const State_ & xBackground,
                                                const State_ & xFirstGuess) {
  Log::trace() << "LinearVariableChange<MODEL>::setTrajectory starting" << std::endl;
  util::Timer timer(classname(), "setTrajectory");
  chvar_->setTrajectory(xBackground.state(), xFirstGuess.state());
  Log::trace() << "LinearVariableChange<MODEL>::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearVariableChange<MODEL>::setTrajectory(const State_ & xBackground) {
  this->setTrajectory(xBackground, xBackground);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearVariableChange<MODEL>::print(std::ostream & os) const {
  Log::trace() << "LinearVariableChange<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *chvar_;
  Log::trace() << "LinearVariableChange<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_LINEARVARIABLECHANGE_H_
