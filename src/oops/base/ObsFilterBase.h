/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSFILTERBASE_H_
#define OOPS_BASE_OBSFILTERBASE_H_

#include <map>
#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Printable.h"

namespace oops {

/// Base class for QC filters applied to observations

// -----------------------------------------------------------------------------

template <typename OBS>
class ObsFilterBase : public util::Printable,
                      private boost::noncopyable {
  typedef GeoVaLs<OBS>             GeoVaLs_;
  typedef ObsDiagnostics<OBS>      ObsDiags_;
  typedef ObsVector<OBS>           ObsVector_;

 public:
  ObsFilterBase() {}
  virtual ~ObsFilterBase() {}

  virtual void preProcess() const = 0;
  virtual void priorFilter(const GeoVaLs_ &) const = 0;
  virtual void postFilter(const ObsVector_ &, const ObsDiags_ &) const = 0;

  virtual Variables requiredVars() const = 0;
  virtual Variables requiredHdiagnostics() const = 0;

 private:
  virtual void print(std::ostream &) const = 0;
};

// =============================================================================

/// ObsFilter Factory
template <typename OBS>
class FilterFactory {
  typedef ObsSpace<OBS>    ObsSpace_;
  template <typename DATA> using ObsDataPtr_ = std::shared_ptr<ObsDataVector<OBS, DATA> >;
 public:
  static std::shared_ptr<ObsFilterBase<OBS>> create(const ObsSpace_ &,
                                                    const eckit::Configuration &,
                                                    ObsDataPtr_<int> flags = ObsDataPtr_<int>(),
                                               ObsDataPtr_<float> obserr = ObsDataPtr_<float>());
  virtual ~FilterFactory() = default;
 protected:
  explicit FilterFactory(const std::string &);
 private:
  virtual ObsFilterBase<OBS> * make(const ObsSpace_ &, const eckit::Configuration &,
                                    ObsDataPtr_<int> &, ObsDataPtr_<float> &) = 0;
  static std::map < std::string, FilterFactory<OBS> * > & getMakers() {
    static std::map < std::string, FilterFactory<OBS> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class OBS, class T>
class FilterMaker : public FilterFactory<OBS> {
  typedef ObsSpace<OBS>    ObsSpace_;
  template <typename DATA> using ObsDataPtr_ = std::shared_ptr<ObsDataVector<OBS, DATA> >;
  virtual ObsFilterBase<OBS> * make(const ObsSpace_ & os, const eckit::Configuration & conf,
                                      ObsDataPtr_<int> & flags, ObsDataPtr_<float> & obserr)
    { return new T(os, conf, flags, obserr); }
 public:
  explicit FilterMaker(const std::string & name) : FilterFactory<OBS>(name) {}
};

// =============================================================================

template <typename OBS>
FilterFactory<OBS>::FilterFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    Log::error() << name << " already registered in obs filter factory." << std::endl;
    ABORT("Element already registered in FilterFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename OBS>
std::shared_ptr<ObsFilterBase<OBS>>
FilterFactory<OBS>::create(const ObsSpace_ & os, const eckit::Configuration & conf,
                             ObsDataPtr_<int> flags, ObsDataPtr_<float> obserr) {
  Log::trace() << "ObsFilterBase<OBS>::create starting" << std::endl;
  const std::string id = conf.getString("filter");
  typename std::map<std::string, FilterFactory<OBS>*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    Log::error() << id << " does not exist in obs filter factory." << std::endl;
    Log::error() << "Obs Filter Factory has " << getMakers().size() << " elements:" << std::endl;
    for (typename std::map<std::string, FilterFactory<OBS>*>::const_iterator
         jj = getMakers().begin(); jj != getMakers().end(); ++jj) {
       Log::error() << "A " << jj->first << " Filter" << std::endl;
    }
    ABORT("Element does not exist in FilterFactory.");
  }
  std::shared_ptr<ObsFilterBase<OBS>> ptr(jloc->second->make(os, conf, flags, obserr));
  Log::trace() << "ObsFilterBase<OBS>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSFILTERBASE_H_
