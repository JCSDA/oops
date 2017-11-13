/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_BASE_FILTERBASE_H_
#define OOPS_BASE_FILTERBASE_H_

#include <boost/noncopyable.hpp>

#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/ObsVector.h"
#include "util/Printable.h"

namespace oops {

/// Base class for QC filters applied to observations

// -----------------------------------------------------------------------------

template <typename MODEL>
class FilterBase : public util::Printable,
                   private boost::noncopyable {
  typedef GeoVaLs<MODEL>             GeoVaLs_;
  typedef ObservationSpace<MODEL>    ObsSpace_;
  typedef ObsVector<MODEL>           ObsVector_;

 public:
  FilterBase() {}
  virtual ~FilterBase() {}

  virtual void postFilter(const GeoVaLs_ &, const ObsVector_ &, const ObsSpace_ &) const =0;

 private:
  virtual void print(std::ostream &) const =0;
};

// =============================================================================

/// ObsFilter Factory
template <typename MODEL>
class FilterFactory {
 public:
  static FilterBase<MODEL> * create(const eckit::Configuration &);
  virtual ~FilterFactory() { getMakers().clear(); }
 protected:
  explicit FilterFactory(const std::string &);
 private:
  virtual FilterBase<MODEL> * make(const eckit::Configuration &) =0;
  static std::map < std::string, FilterFactory<MODEL> * > & getMakers() {
    static std::map < std::string, FilterFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class FilterMaker : public FilterFactory<MODEL> {
  virtual FilterBase<MODEL> * make(const eckit::Configuration & conf)
    { return new T(conf); }
 public:
  explicit FilterMaker(const std::string & name) : FilterFactory<MODEL>(name) {}
};

// =============================================================================

template <typename MODEL>
FilterFactory<MODEL>::FilterFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    Log::error() << name << " already registered in obs filter factory." << std::endl;
    ABORT("Element already registered in FilterFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
FilterBase<MODEL>* FilterFactory<MODEL>::create(const eckit::Configuration & conf) {
  Log::trace() << "FilterBase<MODEL>::create starting" << std::endl;
  Log::debug() << "FilterBase<MODEL>::create conf" << conf << std::endl;
  const std::string id = conf.getString("Filter");
  typename std::map<std::string, FilterFactory<MODEL>*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    Log::error() << id << " does not exist in obs filter factory." << std::endl;
    ABORT("Element does not exist in FilterFactory.");
  }
  FilterBase<MODEL> * ptr = jloc->second->make(conf);
  Log::trace() << "FilterBase<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_FILTERBASE_H_
