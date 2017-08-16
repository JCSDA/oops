/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_INTERFACE_LOCALIZATIONBASE_H_
#define OOPS_INTERFACE_LOCALIZATIONBASE_H_

#include <boost/noncopyable.hpp>
#include <map>
#include <string>

#include "util/Logger.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "eckit/config/Configuration.h"
#include "util/abor1_cpp.h"
#include "util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Base class for localizations

template<typename MODEL>
class LocalizationBase : public util::Printable,
                         private boost::noncopyable {
  typedef Increment<MODEL>        Increment_;
  typedef State<MODEL>            State_;

 public:
  LocalizationBase() {}
  virtual ~LocalizationBase() {}

  void multiply(Increment_ &) const;

 private:
  virtual void multiply(typename MODEL::Increment &) const =0;
  virtual void print(std::ostream &) const =0;
};

// =============================================================================

/// LocalizationFactory Factory
template <typename MODEL>
class LocalizationFactory {
  typedef State<MODEL> State_;
 public:
  static LocalizationBase<MODEL> * create(const State_ &, const eckit::Configuration &);
  virtual ~LocalizationFactory() { getMakers().clear(); }
 protected:
  explicit LocalizationFactory(const std::string &);
 private:
  virtual LocalizationBase<MODEL> * make(const State_ &, const eckit::Configuration &) =0;
  static std::map < std::string, LocalizationFactory<MODEL> * > & getMakers() {
    static std::map < std::string, LocalizationFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class LocalizationMaker : public LocalizationFactory<MODEL> {
  typedef State<MODEL> State_;
  virtual LocalizationBase<MODEL> * make(const State_ & xx, const eckit::Configuration & conf)
    { return new T(xx.state(), conf); }
 public:
  explicit LocalizationMaker(const std::string & name) : LocalizationFactory<MODEL>(name) {}
};

// =============================================================================

template <typename MODEL>
LocalizationFactory<MODEL>::LocalizationFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    Log::error() << name << " already registered in localization factory." << std::endl;
    ABORT("Element already registered in LocalizationFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
LocalizationBase<MODEL>* LocalizationFactory<MODEL>::create(const State_ & xx,
                                                            const eckit::Configuration & conf) {
  Log::trace() << "LocalizationBase<MODEL>::create starting" << std::endl;
  const std::string id = conf.getString("localization");
  typename std::map<std::string, LocalizationFactory<MODEL>*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    Log::error() << id << " does not exist in localization factory." << std::endl;
    ABORT("Element does not exist in LocalizationFactory.");
  }
  LocalizationBase<MODEL> * ptr = jloc->second->make(xx, conf);
  Log::trace() << "LocalizationBase<MODEL>::create done" << std::endl;
  return ptr;
}

// =============================================================================

template <typename MODEL>
void LocalizationBase<MODEL>::multiply(Increment_ & dx) const {
  Log::trace() << "LocalizationBase<MODEL>::multiply starting" << std::endl;
  this->multiply(dx.increment());
  Log::trace() << "LocalizationBase<MODEL>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_LOCALIZATIONBASE_H_
