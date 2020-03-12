/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSLOCALIZATIONBASE_H_
#define OOPS_BASE_OBSLOCALIZATIONBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "oops/interface/ObsSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Base class for generic localizations

template<typename MODEL>
class ObsLocalizationBase : public util::Printable,
                            private boost::noncopyable {
  typedef ObsVector<MODEL> ObsVector_;
 public:
  ObsLocalizationBase() {}
  virtual ~ObsLocalizationBase() {}

  virtual void multiply(ObsVector_ &) const = 0;
};

// =============================================================================

/// ObsLocalizationFactory Factory
template <typename MODEL>
class ObsLocalizationFactory {
  typedef ObsSpace<MODEL>                         ObsSpace_;
 public:
  static std::unique_ptr<ObsLocalizationBase<MODEL>> create(const eckit::Configuration &,
                                                            const ObsSpace_ &);
 protected:
  explicit ObsLocalizationFactory(const std::string &);
 private:
  virtual ObsLocalizationBase<MODEL> * make(const eckit::Configuration &,
                                            const ObsSpace_ &) = 0;
  static std::map < std::string, ObsLocalizationFactory<MODEL> * > & getMakers() {
    static std::map < std::string, ObsLocalizationFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class ObsLocalizationMaker : public ObsLocalizationFactory<MODEL> {
  typedef ObsSpace<MODEL>                         ObsSpace_;
  virtual ObsLocalizationBase<MODEL> * make(const eckit::Configuration & conf,
                                            const ObsSpace_ & obsspace)
    { return new T(conf, obsspace); }
 public:
  explicit ObsLocalizationMaker(const std::string & name) :
    ObsLocalizationFactory<MODEL>(name) {}
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsLocalizationFactory<MODEL>::ObsLocalizationFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    Log::error() << name << " already registered in generic localization factory." << std::endl;
    ABORT("Element already registered in ObsLocalizationFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
std::unique_ptr<ObsLocalizationBase<MODEL>> ObsLocalizationFactory<MODEL>::create(
                              const eckit::Configuration & conf, const ObsSpace_ & obsspace) {
  Log::trace() << "ObsLocalizationBase<MODEL>::create starting" << std::endl;
  const std::string id = conf.getString("localization");
  typename std::map<std::string, ObsLocalizationFactory<MODEL>*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    Log::error() << id << " does not exist in obs localization factory." << std::endl;
    ABORT("Element does not exist in ObsLocalizationFactory.");
  }
  std::unique_ptr<ObsLocalizationBase<MODEL>> ptr(jloc->second->make(conf, obsspace));
  Log::trace() << "ObsLocalizationBase<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSLOCALIZATIONBASE_H_
