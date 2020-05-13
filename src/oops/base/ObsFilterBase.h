/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSFILTERBASE_H_
#define OOPS_BASE_OBSFILTERBASE_H_

#include <map>
#include <string>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

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

template <typename MODEL>
class ObsFilterBase : public util::Printable,
                      private boost::noncopyable {
  typedef GeoVaLs<MODEL>             GeoVaLs_;
  typedef ObsDiagnostics<MODEL>      ObsDiags_;
  typedef ObsVector<MODEL>           ObsVector_;

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
template <typename MODEL>
class FilterFactory {
  typedef ObsSpace<MODEL>    ObsSpace_;
  template <typename DATA> using ObsDataPtr_ = boost::shared_ptr<ObsDataVector<MODEL, DATA> >;
 public:
  static boost::shared_ptr<ObsFilterBase<MODEL>> create(const ObsSpace_ &,
                                                        const eckit::Configuration &,
                                                    ObsDataPtr_<int> flags = ObsDataPtr_<int>(),
                                                 ObsDataPtr_<float> obserr = ObsDataPtr_<float>());
  virtual ~FilterFactory() = default;
 protected:
  explicit FilterFactory(const std::string &);
 private:
  virtual ObsFilterBase<MODEL> * make(const ObsSpace_ &, const eckit::Configuration &,
                                      ObsDataPtr_<int> &, ObsDataPtr_<float> &) = 0;
  static std::map < std::string, FilterFactory<MODEL> * > & getMakers() {
    static std::map < std::string, FilterFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class FilterMaker : public FilterFactory<MODEL> {
  typedef ObsSpace<MODEL>    ObsSpace_;
  template <typename DATA> using ObsDataPtr_ = boost::shared_ptr<ObsDataVector<MODEL, DATA> >;
  virtual ObsFilterBase<MODEL> * make(const ObsSpace_ & os, const eckit::Configuration & conf,
                                      ObsDataPtr_<int> & flags, ObsDataPtr_<float> & obserr)
    { return new T(os, conf, flags, obserr); }
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
boost::shared_ptr<ObsFilterBase<MODEL>>
FilterFactory<MODEL>::create(const ObsSpace_ & os, const eckit::Configuration & conf,
                             ObsDataPtr_<int> flags, ObsDataPtr_<float> obserr) {
  Log::trace() << "ObsFilterBase<MODEL>::create starting" << std::endl;
  const std::string id = conf.getString("Filter");
  typename std::map<std::string, FilterFactory<MODEL>*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    Log::error() << id << " does not exist in obs filter factory." << std::endl;
    Log::error() << "Obs Filter Factory has " << getMakers().size() << " elements:" << std::endl;
    for (typename std::map<std::string, FilterFactory<MODEL>*>::const_iterator
         jj = getMakers().begin(); jj != getMakers().end(); ++jj) {
       Log::error() << "A " << jj->first << " Filter" << std::endl;
    }
    ABORT("Element does not exist in FilterFactory.");
  }
  boost::shared_ptr<ObsFilterBase<MODEL>> ptr(jloc->second->make(os, conf, flags, obserr));
  Log::trace() << "ObsFilterBase<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSFILTERBASE_H_
