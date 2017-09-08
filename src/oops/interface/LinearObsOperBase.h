/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_LINEAROBSOPERBASE_H_
#define OOPS_INTERFACE_LINEAROBSOPERBASE_H_

#include <map>
#include <string>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/interface/LinearObsOperBase.h"
#include "util/abor1_cpp.h"
#include "util/Logger.h"
#include "util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class LinearObsOperBase : public util::Printable,
                          private boost::noncopyable {
  typedef typename MODEL::ModelAtLocations    ModelAtLocations_;
  typedef typename MODEL::ObsAuxControl       ObsAuxControl_;
  typedef typename MODEL::ObsAuxIncrement     ObsAuxIncrement_;
  typedef typename MODEL::ObsVector           ObsVector_;
  typedef typename MODEL::Variables           Variables_;

 public:
  LinearObsOperBase() {}
  virtual ~LinearObsOperBase() {}

/// Obs Operators
  virtual void setTrajectory(const ModelAtLocations_ &, const ObsAuxControl_ &) =0;
  virtual void obsEquivTL(const ModelAtLocations_ &, ObsVector_ &, const ObsAuxIncrement_ &) const =0;
  virtual void obsEquivAD(ModelAtLocations_ &, const ObsVector_ &, ObsAuxIncrement_ &) const =0;

/// Other
  virtual boost::shared_ptr<const Variables_> variables() const =0;  // Required from Model

 private:
  virtual void print(std::ostream &) const =0;
};

// =============================================================================

/// Obs Operator Factory
template <typename MODEL>
class LinearObsOperFactory {
  typedef typename MODEL::ObsSpace ObsSpace_;
 public:
  static LinearObsOperBase<MODEL> * create(const ObsSpace_ &, const eckit::Configuration &);
  virtual ~LinearObsOperFactory() { getMakers().clear(); }
 protected:
  explicit LinearObsOperFactory(const std::string &);
 private:
  virtual LinearObsOperBase<MODEL> * make(const ObsSpace_ &, const eckit::Configuration &) =0;
  static std::map < std::string, LinearObsOperFactory<MODEL> * > & getMakers() {
    static std::map < std::string, LinearObsOperFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class LinearObsOpMaker : public LinearObsOperFactory<MODEL> {
  typedef typename MODEL::ObsSpace const ObsSpace_;
  virtual LinearObsOperBase<MODEL> * make(const ObsSpace_ & odb,
                                          const eckit::Configuration & conf)
    { return new T(odb, conf); }
 public:
  explicit LinearObsOpMaker(const std::string & name) : LinearObsOperFactory<MODEL>(name) {}
};

// =============================================================================

template <typename MODEL>
LinearObsOperFactory<MODEL>::LinearObsOperFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    Log::error() << name << " already registered in LinearObsOperFactory." << std::endl;
    ABORT("Element already registered in LinearObsOperFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearObsOperBase<MODEL>* LinearObsOperFactory<MODEL>::create(const ObsSpace_ & odb,
                                                              const eckit::Configuration & conf) {
  const std::string id = conf.getString("ObsType");
  typename std::map<std::string, LinearObsOperFactory<MODEL>*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    Log::error() << id << " does not exist in LinearObsOperFactory." << std::endl;
    ABORT("Element does not exist in LinearObsOperFactory.");
  }
  LinearObsOperBase<MODEL> * ptr = jloc->second->make(odb, conf);
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_LINEAROBSOPERBASE_H_
