/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_VARIABLECHANGEBASE_H_
#define OOPS_BASE_VARIABLECHANGEBASE_H_

#include <map>
#include <memory>
#include <string>
#include <boost/noncopyable.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/State.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

/// Base class for generic variable transform

template <typename MODEL>
class VariableChangeBase : public util::Printable,
                           private boost::noncopyable {
  typedef State<MODEL>               State_;

 public:
  explicit VariableChangeBase(const eckit::Configuration &);
  virtual ~VariableChangeBase() {}

  void setInputVariables(const Variables & vars) { varin_.reset(new Variables(vars)); }
  void setOutputVariables(const Variables & vars) { varout_.reset(new Variables(vars)); }

  virtual void changeVar(const State_ &, State_ &) const = 0;
  virtual void changeVarInverse(const State_ &, State_ &) const = 0;

  State_ changeVar(const State_ &) const;
  State_ changeVarInverse(const State_ &) const;

 private:
  virtual void print(std::ostream &) const = 0;
  std::unique_ptr<Variables> varin_;
  std::unique_ptr<Variables> varout_;
};

// -----------------------------------------------------------------------------

/// VariableChangeFactory Factory
template <typename MODEL>
class VariableChangeFactory {
  typedef Geometry<MODEL>   Geometry_;
 public:
  static VariableChangeBase<MODEL> * create(const eckit::Configuration &, const Geometry_ &);
  virtual ~VariableChangeFactory() { getMakers().clear(); }
 protected:
  explicit VariableChangeFactory(const std::string &);
 private:
  virtual VariableChangeBase<MODEL> * make(const eckit::Configuration &, const Geometry_ &) = 0;
  static std::map < std::string, VariableChangeFactory<MODEL> * > & getMakers() {
    static std::map < std::string, VariableChangeFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class VariableChangeMaker : public VariableChangeFactory<MODEL> {
  typedef Geometry<MODEL>   Geometry_;
  virtual VariableChangeBase<MODEL> * make(const eckit::Configuration & conf,
                                           const Geometry_ & resol)
    { return new T(resol, conf); }
 public:
  explicit VariableChangeMaker(const std::string & name)
    : VariableChangeFactory<MODEL>(name) {}
};

// -----------------------------------------------------------------------------

template <typename MODEL>
VariableChangeFactory<MODEL>::VariableChangeFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    Log::error() << name << " already registered in the variable change factory."  << std::endl;
    ABORT("Element already registered in VariableChangeFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
VariableChangeBase<MODEL> * VariableChangeFactory<MODEL>::create(const eckit::Configuration & conf,
                                                                 const Geometry_ & resol)
{
  Log::trace() << "VariableChangeBase<MODEL>::create starting" << std::endl;
  std::string id;
  if (conf.has("varchange")) {
    id = conf.getString("varchange");
  } else {
    id = "Identity";
  }
  typename std::map<std::string, VariableChangeFactory<MODEL>*>::iterator
    jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    Log::error() << id << " does not exist in the variable change factory factory." << std::endl;
    ABORT("Element does not exist in VariableChangeFactory.");
  }
  VariableChangeBase<MODEL> * ptr = jerr->second->make(conf, resol);
  Log::trace() << "VariableChangeBase<MODEL>::create done" << std::endl;
  return ptr;
}

// =============================================================================

template<typename MODEL>
VariableChangeBase<MODEL>::VariableChangeBase(const eckit::Configuration & conf)
  : varin_(), varout_()
{
  if (conf.has("inputVariables")) {
    varin_.reset(new Variables(conf.getSubConfiguration("inputVariables")));
    Log::trace() << "VariableChangeBase<MODEL>::VariableChangeBase inputvars: "
                 << *varin_ << std::endl;
  }
  if (conf.has("outputVariables")) {
    varout_.reset(new Variables(conf.getSubConfiguration("outputVariables")));
    Log::trace() << "VariableChangeBase<MODEL>::VariableChangeBase outputvars: "
                 << *varout_ << std::endl;
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL> VariableChangeBase<MODEL>::changeVar(const State_ & dxin) const {
  ASSERT(varin_);
  State_ dxout(dxin.geometry(), *varin_, dxin.validTime());
  this->multiply(dxin, dxout);
  return dxout;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL> VariableChangeBase<MODEL>::changeVarInverse(const State_ & dxin) const {
  ASSERT(varout_);
  State_ dxout(dxin.geometry(), *varout_, dxin.validTime());
  this->multiplyInverse(dxin, dxout);
  return dxout;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_VARIABLECHANGEBASE_H_
