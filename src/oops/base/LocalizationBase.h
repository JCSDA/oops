/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_LOCALIZATIONBASE_H_
#define OOPS_BASE_LOCALIZATIONBASE_H_

#include <map>
#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "oops/assimilation/Increment4D.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// Model-space localization base class
template <typename MODEL>
class LocalizationBase : public util::Printable,
                         private boost::noncopyable {
  typedef Increment<MODEL>                        Increment_;
  typedef Increment4D<MODEL>                      Increment4D_;

 public:
  LocalizationBase() = default;
  virtual ~LocalizationBase() = default;

  virtual void multiply(Increment_ &) const = 0;
  /// default implementation of 4D localization (used in model-specific implementations)
  virtual void multiply(Increment4D_ &) const;

 private:
  virtual void print(std::ostream &) const = 0;
};

// -----------------------------------------------------------------------------
/// Default implementation of 4D localization (used in model-specific implementations)
template <typename MODEL>
void LocalizationBase<MODEL>::multiply(Increment4D_ & dx) const {
  Log::trace() << "LocalizationBase<MODEL>::multiply starting" << std::endl;
  // Sum over timeslots
  Increment_ dxtmp(dx[dx.first()]);
  for (int isub = dx.first()+1; isub <= dx.last(); ++isub) {
     dxtmp.axpy(1.0, dx[isub], false);
  }

  // Apply 3D localization
  this->multiply(dxtmp);

  // Copy result to all timeslots
  for (int isub = dx.first(); isub <= dx.last(); ++isub) {
     dx[isub].zero();
     dx[isub].axpy(1.0, dxtmp, false);
  }
  Log::trace() << "LocalizationBase<MODEL>::multiply done" << std::endl;
}


// =============================================================================

/// Localization Factory
template <typename MODEL>
class LocalizationFactory {
  typedef Geometry<MODEL>                             Geometry_;
 public:
  static std::unique_ptr<LocalizationBase<MODEL>> create(const Geometry_ &,
                          const eckit::Configuration &);
  virtual ~LocalizationFactory() = default;
 protected:
  explicit LocalizationFactory(const std::string &);
 private:
  virtual LocalizationBase<MODEL> * make(const Geometry_ &,
                                         const eckit::Configuration &) = 0;
  static std::map < std::string, LocalizationFactory<MODEL> * > & getMakers() {
    static std::map < std::string, LocalizationFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class LocalizationMaker : public LocalizationFactory<MODEL> {
  typedef Geometry<MODEL>                             Geometry_;
  virtual LocalizationBase<MODEL> * make(const Geometry_ & geometry,
                                         const eckit::Configuration & conf)
    { return new T(geometry, conf); }
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
std::unique_ptr<LocalizationBase<MODEL>>
LocalizationFactory<MODEL>::create(const Geometry_ & geometry,
                                   const eckit::Configuration & conf) {
  Log::trace() << "LocalizationBase<MODEL>::create starting" << std::endl;
  const std::string id = conf.getString("localization method");
  typename std::map<std::string, LocalizationFactory<MODEL>*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    Log::error() << id << " does not exist in localization factory." << std::endl;
    Log::error() << "Obs Localization Factory has " << getMakers().size()
                 << " elements:" << std::endl;
    for (typename std::map<std::string, LocalizationFactory<MODEL>*>::const_iterator
         jj = getMakers().begin(); jj != getMakers().end(); ++jj) {
       Log::error() << "A " << jj->first << " Localization" << std::endl;
    }
    ABORT("Element does not exist in LocalizationFactory.");
  }
  std::unique_ptr<LocalizationBase<MODEL>> ptr(jloc->second->make(geometry, conf));
  Log::trace() << "LocalizationBase<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_LOCALIZATIONBASE_H_
