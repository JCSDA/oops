/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_STAT_VARIABLECHANGE_H_
#define OOPS_STAT_VARIABLECHANGE_H_

#include <boost/noncopyable.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

/// Derived class of generic variable transform for 

template <typename MODEL>
class StatVariableChange : public VariableChangeBase<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

public:
  static const std::string classname() {return "oops::StatVariableChange";}

  explicit StatVariableChange(const eckit::Configuration &);
  virtual ~StatVariableChange();

  void linearize(const State_ &,  const Geometry_ &) override;
  void transform(const Increment_ &, Increment_ &) const override;
  void transformInverse(const Increment_ &, Increment_ &) const override;
  void transformAD(const Increment_ &, Increment_ &) const override;
  void transformInverseAD(const Increment_ &, Increment_ &) const override;

private:
  void print(std::ostream &) const override;


  //StatVarData populate(const varin_ , const varout_, const eckit::Configuration);

  
};

template<typename MODEL>
StatVariableChange<MODEL>::StatVariableChange(const eckit::Configuration & conf)
  : VariableChangeBase<MODEL>(conf)
{
  Log::trace() << "StatVariableChange<MODEL>::StatVariableChange starting" << std::endl;
  Log::trace() << "StatVariableChange<MODEL>::StatVariableChange done" << std::endl;
}

template<typename MODEL>
StatVariableChange<MODEL>::~StatVariableChange() {
  Log::trace() << "StatVariableChange<MODEL>::~StatVariableChange starting" << std::endl;
  Log::trace() << "StatVariableChange<MODEL>::~StatVariableChange done" << std::endl;
}

template<typename MODEL>
void StatVariableChange<MODEL>::linearize(const State_ & x,  const Geometry_ & resol) {
  Log::trace() << "StatVariableChange<MODEL>::linearize starting" << std::endl;
  Log::trace() << "StatVariableChange<MODEL>::linearize done" << std::endl;
}

template<typename MODEL>
void StatVariableChange<MODEL>::transform(const Increment_ & in, Increment_ & out) const {
  Log::trace() << "StatVariableChange<MODEL>::transform starting" << std::endl;

  UnstructuredGrid ug;
  in.convert_to(ug);
  out.convert_from(ug);

  Log::trace() << "StatVariableChange<MODEL>::transform done" << std::endl;
}

template<typename MODEL>
void StatVariableChange<MODEL>::transformInverse(const Increment_ & in, Increment_ & out) const {
  Log::trace() << "StatVariableChange<MODEL>::transformInverse starting" << std::endl;

  UnstructuredGrid ug;
  in.convert_to(ug);
  out.convert_from(ug);

  Log::trace() << "StatVariableChange<MODEL>::transformInverse done" << std::endl;
}

template<typename MODEL>
void StatVariableChange<MODEL>::transformAD(const Increment_ & in, Increment_ & out) const {
  Log::trace() << "StatVariableChange<MODEL>::transformAD starting" << std::endl;

  UnstructuredGrid ug;
  in.convert_to(ug);
  out.convert_from(ug);

  Log::trace() << "StatVariableChange<MODEL>::transformAD done" << std::endl;
}

template<typename MODEL>
void StatVariableChange<MODEL>::transformInverseAD(const Increment_ & in, Increment_ & out) const {
  Log::trace() << "StatVariableChange<MODEL>::transformInverseAD starting" << std::endl;

  UnstructuredGrid ug;
  in.convert_to(ug);
  out.convert_from(ug);

  Log::trace() << "StatVariableChange<MODEL>::transformInverseAD done" << std::endl;
}

template<typename MODEL>
void StatVariableChange<MODEL>::print(std::ostream & os) const {
}

}  // namespace oops

#endif  // OOPS_STAT_VARIABLECHANGE_H_
