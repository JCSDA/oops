/*
 * (C) Copyright 2018-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_CHANGEVARIABLES_H_
#define OOPS_INTERFACE_CHANGEVARIABLES_H_

#include <memory>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/State.h"
#include "oops/interface/VariableChange.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

/// \brief Encapsulates the nonlinear variable change
/// There should not be a factory for ChangeVariable, it should be a trait class.
/// This is a temporary solution to get around it until we can fix the models.
template <typename MODEL>
class ChangeVariables : public util::Printable,
                        private util::ObjectCounter<ChangeVariables<MODEL> >  {
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;
  typedef VariableChange<MODEL>      VariableChange_;

 public:
  static const std::string classname() {return "oops::ChangeVariables";}

  ChangeVariables(const eckit::Configuration &, const Geometry_ &,
                  const Variables &, const Variables &);
  virtual ~ChangeVariables();
  ChangeVariables(const ChangeVariables &) = delete;
  ChangeVariables(ChangeVariables &&) = default;
  const ChangeVariables & operator=(const ChangeVariables &) = delete;
  ChangeVariables & operator=(ChangeVariables &&) = default;

  /// change variable from state \p xin to \p xout
  void changeVar(const State_ & xin, State_ & xout) const;
  /// inverse of changeVar, change variables back from \p xout to \p xin
  void changeVarInverse(const State_ & xout, State_ & xin) const;

 private:
  void print(std::ostream &) const override;

  std::unique_ptr<VariableChange_> chvar_;  /// pointer to the VariableChange implementation
};

// =============================================================================

template<typename MODEL>
ChangeVariables<MODEL>::ChangeVariables(const eckit::Configuration & conf, const Geometry_ & geom,
                                        const Variables & varin, const Variables & varout)
  : chvar_()
{
  Log::trace() << "ChangeVariables<MODEL>::ChangeVariables starting" << std::endl;
//  Comment timer for now to avoid double counting
//  util::Timer timer(classname(), "ChangeVariables");
  eckit::LocalConfiguration chconf(conf);
  if (!chconf.has("variable change")) chconf.set("variable change", "default");
  chconf.set("input variables", varin.variables());
  chconf.set("output variables", varout.variables());
  chvar_.reset(new VariableChange_(geom, chconf));
  Log::trace() << "ChangeVariables<MODEL>::ChangeVariables done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ChangeVariables<MODEL>::~ChangeVariables() {
  Log::trace() << "ChangeVariables<MODEL>::~ChangeVariables starting" << std::endl;
//  util::Timer timer(classname(), "~ChangeVariables");
  chvar_.reset();
  Log::trace() << "ChangeVariables<MODEL>::~ChangeVariables done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ChangeVariables<MODEL>::changeVar(const State_ & x1, State_ & x2) const {
  Log::trace() << "ChangeVariables<MODEL>::changeVar starting" << std::endl;
//  util::Timer timer(classname(), "changeVar");
  chvar_->changeVar(x1, x2);
  Log::trace() << "ChangeVariables<MODEL>::changeVar done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ChangeVariables<MODEL>::changeVarInverse(const State_ & x1, State_ & x2) const {
  Log::trace() << "ChangeVariables<MODEL>::changeVarInverse starting" << std::endl;
//  util::Timer timer(classname(), "changeVarInverse");
  chvar_->changeVarInverse(x1, x2);
  Log::trace() << "ChangeVariables<MODEL>::changeVarInverse done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ChangeVariables<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ChangeVariables<MODEL>::print starting" << std::endl;
//  util::Timer timer(classname(), "print");
  os << *chvar_;
  Log::trace() << "ChangeVariables<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_CHANGEVARIABLES_H_
