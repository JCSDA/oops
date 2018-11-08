/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_LINEARVARIABLECHANGEBASE_H_
#define OOPS_BASE_LINEARVARIABLECHANGEBASE_H_

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
class LinearVariableChangeBase : public util::Printable,
                                 private boost::noncopyable {
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
  explicit LinearVariableChangeBase(const eckit::Configuration &);
  virtual ~LinearVariableChangeBase() {}

  void setInputVariables(const Variables & vars) { varin_.reset(new Variables(vars)); }
  void setOutputVariables(const Variables & vars) { varout_.reset(new Variables(vars)); }

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

/// LinearVariableChangeFactory Factory
template <typename MODEL>
class LinearVariableChangeFactory {
  typedef Geometry<MODEL>   Geometry_;
  typedef State<MODEL>      State_;
 public:
  static LinearVariableChangeBase<MODEL> * create(const State_ &, const State_ &,
                                                  const Geometry_ &,
                                                  const eckit::Configuration &);
  virtual ~LinearVariableChangeFactory() { getMakers().clear(); }
 protected:
  explicit LinearVariableChangeFactory(const std::string &);
 private:
  virtual LinearVariableChangeBase<MODEL> * make(const State_ &, const State_ &,
                                                 const Geometry_ &,
                                                 const eckit::Configuration &) = 0;
  static std::map < std::string, LinearVariableChangeFactory<MODEL> * > & getMakers() {
    static std::map < std::string, LinearVariableChangeFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class LinearVariableChangeMaker : public LinearVariableChangeFactory<MODEL> {
  typedef Geometry<MODEL>   Geometry_;
  typedef State<MODEL>      State_;
  virtual LinearVariableChangeBase<MODEL> * make(const State_ & bg, const State_ & fg,
                                                 const Geometry_ & geom,
                                                 const eckit::Configuration & conf)
    { return new T(bg, fg, geom, conf); }
 public:
  explicit LinearVariableChangeMaker(const std::string & name)
    : LinearVariableChangeFactory<MODEL>(name) {}
};

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearVariableChangeFactory<MODEL>::LinearVariableChangeFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    Log::error() << name << " already registered in the variable change factory."  << std::endl;
    ABORT("Element already registered in LinearVariableChangeFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearVariableChangeBase<MODEL> * LinearVariableChangeFactory<MODEL>::create(
     const State_ & bg, const State_ & fg,
     const Geometry_ & geom, const eckit::Configuration & conf) {
  Log::trace() << "LinearVariableChangeBase<MODEL>::create starting" << std::endl;
  const std::string id = conf.getString("varchange");
  typename std::map<std::string, LinearVariableChangeFactory<MODEL>*>::iterator
    jerr = getMakers().find(id);

  if (jerr == getMakers().end()) {
    Log::error() << id << " does not exist in LinearVariableChangeFactory." << std::endl;
    Log::error() << "Factory contains " << getMakers().size() << " elements:" << std::endl;
    for (typename std::map<std::string, LinearVariableChangeFactory<MODEL>*>::const_iterator
         jj = getMakers().begin(); jj !=  getMakers().end(); ++jj) {
      Log::error() << "A " << jj->first << " variable change option" << std::endl;
    }
    ABORT("Element does not exist in LinearVariableChangeFactory.");
  }

  LinearVariableChangeBase<MODEL> * ptr = jerr->second->make(bg, fg, geom, conf);
  Log::trace() << "LinearVariableChangeBase<MODEL>::create done" << std::endl;
  return ptr;
}

// =============================================================================

template<typename MODEL>
LinearVariableChangeBase<MODEL>::LinearVariableChangeBase(const eckit::Configuration & conf)
  : varin_(), varout_()
{
  if (conf.has("inputVariables")) {
    varin_.reset(new Variables(conf.getSubConfiguration("inputVariables")));
    Log::trace() << "LinearVariableChangeBase<MODEL>::LinearVariableChangeBase inputvars: "
                 << *varin_ << std::endl;
  }
  if (conf.has("outputVariables")) {
    varout_.reset(new Variables(conf.getSubConfiguration("outputVariables")));
    Log::trace() << "LinearVariableChangeBase<MODEL>::LinearVariableChangeBase outputvars: "
                 << *varout_ << std::endl;
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> LinearVariableChangeBase<MODEL>::multiply(const Increment<MODEL> & dxin) const {
  ASSERT(varin_);
  Increment_ dxout(dxin.geometry(), *varin_, dxin.validTime());
  this->multiply(dxin, dxout);
  return dxout;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> LinearVariableChangeBase<MODEL>::multiplyAD(const Increment_ & dxin) const {
  ASSERT(varout_);
  Increment_ dxout(dxin.geometry(), *varout_, dxin.validTime());
  this->multiplyAD(dxin, dxout);
  return dxout;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> LinearVariableChangeBase<MODEL>::multiplyInverse(const Increment_ & dxin) const {
  ASSERT(varout_);
  Increment_ dxout(dxin.geometry(), *varout_, dxin.validTime());
  this->multiplyInverse(dxin, dxout);
  return dxout;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> LinearVariableChangeBase<MODEL>::multiplyInverseAD(const Increment_ & dxin) const {
  ASSERT(varin_);
  Increment_ dxout(dxin.geometry(), *varin_, dxin.validTime());
  this->multiplyInverseAD(dxin, dxout);
  return dxout;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_LINEARVARIABLECHANGEBASE_H_
