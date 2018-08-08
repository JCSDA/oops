/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_CHANGEVARIABLE_H_
#define OOPS_INTERFACE_CHANGEVARIABLE_H_

#include <string>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/generic/VariableChangeBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// Wrapper for change of variable

template <typename MODEL>
class ChangeVariable : public oops::VariableChangeBase<MODEL>,
                       public util::Printable {
  typedef typename MODEL::ChangeVar  ChangeVar_;
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
  static const std::string classname() {return "oops::ChangeVariable";}

  explicit ChangeVariable(const eckit::Configuration &);
  virtual ~ChangeVariable();

  void linearize(const State_ &, const Geometry_ &) override;

  void transform(const Increment_ &, Increment_ &) const override;
  void transformInverse(const Increment_ &, Increment_ &) const override;
  void transformAD(const Increment_ &, Increment_ &) const override;
  void transformInverseAD(const Increment_ &, Increment_ &) const override;
 
 private:
  void print(std::ostream &) const override;

  boost::scoped_ptr<ChangeVar_> chvar_;
};

// =============================================================================

template<typename MODEL>
ChangeVariable<MODEL>::ChangeVariable(const eckit::Configuration & conf)
  : VariableChangeBase<MODEL>(conf), chvar_()
{
  Log::trace() << "ChangeVariable<MODEL>::ChangeVariable starting" << std::endl;
  util::Timer timer(classname(), "ChangeVariable");
  chvar_.reset(new ChangeVar_(conf));
  Log::trace() << "ChangeVariable<MODEL>::ChangeVariable done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ChangeVariable<MODEL>::~ChangeVariable() {
  Log::trace() << "ChangeVariable<MODEL>::~ChangeVariable starting" << std::endl;
  util::Timer timer(classname(), "~ChangeVariable");
  chvar_.reset();
  Log::trace() << "ChangeVariable<MODEL>::~ChangeVariable done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ChangeVariable<MODEL>::linearize(const State_ & xx, const Geometry_ & resol) {
  Log::trace() << "ChangeVariable<MODEL>::linearize starting" << std::endl;
  util::Timer timer(classname(), "linearize");
  chvar_->linearize(xx.state(), resol.geometry());
  Log::trace() << "ChangeVariable<MODEL>::linearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ChangeVariable<MODEL>::transform(const Increment_ & dx1, Increment_ & dx2) const {
  Log::trace() << "ChangeVariable<MODEL>::transform starting" << std::endl;
  util::Timer timer(classname(), "transform");
  chvar_->transform(dx1.increment(), dx2.increment());
  Log::trace() << "ChangeVariable<MODEL>::transform done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ChangeVariable<MODEL>::transformInverse(const Increment_ & dx1, Increment_ & dx2) const {
  Log::trace() << "ChangeVariable<MODEL>::transformInverse starting" << std::endl;
  util::Timer timer(classname(), "transformInverse");
  chvar_->transformInverse(dx1.increment(), dx2.increment());
  Log::trace() << "ChangeVariable<MODEL>::transformInverse done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ChangeVariable<MODEL>::transformAD(const Increment_ & dx1, Increment_ & dx2) const {
  Log::trace() << "ChangeVariable<MODEL>::transformAD starting" << std::endl;
  util::Timer timer(classname(), "transformAD");
  chvar_->transformAD(dx1.increment(), dx2.increment());
  Log::trace() << "ChangeVariable<MODEL>::transformAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ChangeVariable<MODEL>::transformInverseAD(const Increment_ & dx1, Increment_ & dx2) const {
  Log::trace() << "ChangeVariable<MODEL>::transformInverseAD starting" << std::endl;
  util::Timer timer(classname(), "transformInverseAD");
  chvar_->transformInverseAD(dx1.increment(), dx2.increment());
  Log::trace() << "ChangeVariable<MODEL>::transformInverseAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ChangeVariable<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ChangeVariable<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *chvar_;
  Log::trace() << "ChangeVariable<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_CHANGEVARIABLE_H_
