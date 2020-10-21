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
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Base class for generic localizations

template<typename OBS>
class ObsLocalizationBase : public util::Printable,
                            private boost::noncopyable {
  typedef ObsVector<OBS> ObsVector_;
 public:
  ObsLocalizationBase() {}
  virtual ~ObsLocalizationBase() {}

  virtual void multiply(ObsVector_ &) const = 0;
};

// =============================================================================

/// ObsLocalizationFactory Factory
template <typename OBS>
class ObsLocalizationFactory {
  typedef ObsSpace<OBS>                         ObsSpace_;
 public:
  static std::unique_ptr<ObsLocalizationBase<OBS>> create(const eckit::Configuration &,
                                                            const ObsSpace_ &);
 protected:
  explicit ObsLocalizationFactory(const std::string &);
 private:
  virtual ObsLocalizationBase<OBS> * make(const eckit::Configuration &,
                                            const ObsSpace_ &) = 0;
  static std::map < std::string, ObsLocalizationFactory<OBS> * > & getMakers() {
    static std::map < std::string, ObsLocalizationFactory<OBS> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class OBS, class T>
class ObsLocalizationMaker : public ObsLocalizationFactory<OBS> {
  typedef ObsSpace<OBS>                         ObsSpace_;
  virtual ObsLocalizationBase<OBS> * make(const eckit::Configuration & conf,
                                            const ObsSpace_ & obsspace)
    { return new T(conf, obsspace); }
 public:
  explicit ObsLocalizationMaker(const std::string & name) :
    ObsLocalizationFactory<OBS>(name) {}
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsLocalizationFactory<OBS>::ObsLocalizationFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in obs localization factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename OBS>
std::unique_ptr<ObsLocalizationBase<OBS>> ObsLocalizationFactory<OBS>::create(
                              const eckit::Configuration & conf, const ObsSpace_ & obsspace) {
  Log::trace() << "ObsLocalizationBase<OBS>::create starting" << std::endl;
  const std::string id = conf.getString("localization method");
  typename std::map<std::string, ObsLocalizationFactory<OBS>*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    throw std::runtime_error(id + " does not exist in obs localization factory.");
  }
  std::unique_ptr<ObsLocalizationBase<OBS>> ptr(jloc->second->make(conf, obsspace));
  Log::trace() << "ObsLocalizationBase<OBS>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSLOCALIZATIONBASE_H_
