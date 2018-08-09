/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_VARIABLECHANGEBASE_H_
#define OOPS_GENERIC_VARIABLECHANGEBASE_H_

#include <map>
#include <string>
#include <boost/noncopyable.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
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
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
  explicit VariableChangeBase(const eckit::Configuration &);
  virtual ~VariableChangeBase() {}

  void setInputVariables(Variables & vars) { varin_.reset(new Variables(vars)); }
  void setOutputVariables(Variables & vars) { varout_.reset(new Variables(vars)); }

  virtual void linearize(const State_ &, const Geometry_ &) = 0;

  virtual void multiply(const Increment_ &, Increment_ &) const = 0;
  virtual void multiplyInverse(const Increment_ &, Increment_ &) const = 0;
  virtual void multiplyAD(const Increment_ &, Increment_ &) const = 0;
  virtual void multiplyInverseAD(const Increment_ &, Increment_ &) const = 0;

  Increment_ multiply(const Increment_ &) const;
  Increment_ multiplyInverse(const Increment_ &) const;
  Increment_ multiplyAD(const Increment_ &) const;
  Increment_ multiplyInverseAD(const Increment_ &) const;

 private:
  virtual void print(std::ostream &) const = 0;
  boost::scoped_ptr<Variables> varin_;
  boost::scoped_ptr<Variables> varout_;
};

// -----------------------------------------------------------------------------

/// VariableChangeFactory Factory
template <typename MODEL>
class VariableChangeFactory {
  typedef Geometry<MODEL>   Geometry_;
  typedef State<MODEL>      State_;
 public:
  static VariableChangeBase<MODEL> * create(const eckit::Configuration &);
  virtual ~VariableChangeFactory() { getMakers().clear(); }
 protected:
  explicit VariableChangeFactory(const std::string &);
 private:
  virtual VariableChangeBase<MODEL> * make(const eckit::Configuration &) = 0;
  static std::map < std::string, VariableChangeFactory<MODEL> * > & getMakers() {
    static std::map < std::string, VariableChangeFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class VariableChangeMaker : public VariableChangeFactory<MODEL> {
  typedef Geometry<MODEL>   Geometry_;
  virtual VariableChangeBase<MODEL> * make(const eckit::Configuration & conf)
    { return new T(conf); }
 public:
  explicit VariableChangeMaker(const std::string & name) : VariableChangeFactory<MODEL>(name) {}
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
VariableChangeBase<MODEL>* VariableChangeFactory<MODEL>::create(const eckit::Configuration & conf) {
  Log::trace() << "VariableChangeBase<MODEL>::create starting" << std::endl;
  const std::string id = conf.getString("varchange");
  typename std::map<std::string, VariableChangeFactory<MODEL>*>::iterator
    jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    Log::error() << id << " does not exist in the variable change factory factory." << std::endl;
    ABORT("Element does not exist in VariableChangeFactory.");
  }
  VariableChangeBase<MODEL> * ptr = jerr->second->make(conf);
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
    Log::trace() << "VariableChangeBase<MODEL>::VariableChangeBase inputvars: " << *varin_ << std::endl;
  }
  if (conf.has("outputVariables")) {
    varout_.reset(new Variables(conf.getSubConfiguration("outputVariables")));
    Log::trace() << "VariableChangeBase<MODEL>::VariableChangeBase outputvars: " << *varout_ << std::endl;
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> VariableChangeBase<MODEL>::multiply(const Increment<MODEL> & dxin) const {
  ASSERT(varout_);
  Increment_ dxout(dxin.geometry(), *varout_, dxin.validTime());
  this->multiply(dxin, dxout);
  return dxout;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> VariableChangeBase<MODEL>::multiplyAD(const Increment_ & dxin) const {
  ASSERT(varin_);
  Increment_ dxout(dxin.geometry(), *varin_, dxin.validTime());
  this->multiplyAD(dxin, dxout);
  return dxout;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> VariableChangeBase<MODEL>::multiplyInverse(const Increment_ & dxin) const {
  ASSERT(varin_);
  Increment_ dxout(dxin.geometry(), *varin_, dxin.validTime());
  this->multiplyInverse(dxin, dxout);
  return dxout;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> VariableChangeBase<MODEL>::multiplyInverseAD(const Increment_ & dxin) const {
  ASSERT(varout_);
  Increment_ dxout(dxin.geometry(), *varout_, dxin.validTime());
  this->multiplyInverseAD(dxin, dxout);
  return dxout;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_VARIABLECHANGEBASE_H_
