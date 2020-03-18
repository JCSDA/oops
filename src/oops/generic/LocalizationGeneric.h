/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_LOCALIZATIONGENERIC_H_
#define OOPS_GENERIC_LOCALIZATIONGENERIC_H_

#include <map>
#include <string>
#include <vector>
#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "oops/assimilation/Increment4D.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Base class for generic localizations

template<typename MODEL>
class LocalizationGeneric : public util::Printable,
                            private boost::noncopyable {
  typedef Increment<MODEL>                        Increment_;
  typedef Increment4D<MODEL>                      Increment4D_;

 public:
  LocalizationGeneric() {}
  virtual ~LocalizationGeneric() {}

  void doMultiply(Increment_ &) const;
  void doMultiply(Increment4D_ &) const;

 private:
  virtual void multiply(Increment_ &) const = 0;
  virtual void multiply(Increment4D_ &) const = 0;
  virtual void print(std::ostream &) const = 0;
};

// =============================================================================

/// LocalizationGenericFactory Factory
template <typename MODEL>
class LocalizationGenericFactory {
  typedef Geometry<MODEL>                         Geometry_;
  typedef boost::shared_ptr<IncrementEnsemble<MODEL>> EnsemblePtr_;
 public:
  static LocalizationGeneric<MODEL> * create(const Geometry_ &,
                                             const EnsemblePtr_,
                                             const eckit::Configuration &);
  virtual ~LocalizationGenericFactory() = default;
 protected:
  explicit LocalizationGenericFactory(const std::string &);
 private:
  virtual LocalizationGeneric<MODEL> * make(const Geometry_ &,
                                            const EnsemblePtr_,
                                            const eckit::Configuration &) = 0;
  static std::map < std::string, LocalizationGenericFactory<MODEL> * > & getMakers() {
    static std::map < std::string, LocalizationGenericFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class LocalizationGenericMaker : public LocalizationGenericFactory<MODEL> {
  typedef Geometry<MODEL>                         Geometry_;
  typedef boost::shared_ptr<IncrementEnsemble<MODEL>> EnsemblePtr_;
  virtual LocalizationGeneric<MODEL> * make(const Geometry_ & resol,
                                            const EnsemblePtr_ ens,
                                            const eckit::Configuration & conf)
    { return new T(resol, ens, conf); }
 public:
  explicit LocalizationGenericMaker(const std::string & name) :
    LocalizationGenericFactory<MODEL>(name) {}
};

// -----------------------------------------------------------------------------

template <typename MODEL>
LocalizationGenericFactory<MODEL>::LocalizationGenericFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    Log::error() << name << " already registered in generic localization factory." << std::endl;
    ABORT("Element already registered in LocalizationGenericFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
LocalizationGeneric<MODEL>* LocalizationGenericFactory<MODEL>::create(const Geometry_ & resol,
                                                              const EnsemblePtr_ ens,
                                                              const eckit::Configuration & conf) {
  Log::trace() << "LocalizationGeneric<MODEL>::create starting" << std::endl;
  const std::string id = conf.getString("localization");
  typename std::map<std::string, LocalizationGenericFactory<MODEL>*>::iterator
    jloc = getMakers().find(id);
  LocalizationGeneric<MODEL> * ptr = NULL;
  if (jloc != getMakers().end()) {
     ptr = jloc->second->make(resol, ens, conf);
  }
  Log::trace() << "LocalizationGeneric<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LocalizationGeneric<MODEL>::doMultiply(Increment_ & dx) const {
  Log::trace() << "LocalizationGeneric<MODEL>::multiply starting" << std::endl;
  this->multiply(dx);
  Log::trace() << "LocalizationGeneric<MODEL>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LocalizationGeneric<MODEL>::doMultiply(Increment4D_ & dx) const {
  Log::trace() << "LocalizationGeneric<MODEL>::multiply starting" << std::endl;
  this->multiply(dx);
  Log::trace() << "LocalizationGeneric<MODEL>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_LOCALIZATIONGENERIC_H_
