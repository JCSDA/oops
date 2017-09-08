/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_OBSOPERATORBASE_H_
#define OOPS_INTERFACE_OBSOPERATORBASE_H_

#include <map>
#include <string>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "util/abor1_cpp.h"
#include "util/Logger.h"
#include "util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Base class for observation operators

template <typename MODEL>
class ObsOperatorBase : public util::Printable,
                        private boost::noncopyable {
  typedef typename MODEL::ModelAtLocations    ModelAtLocations_;
  typedef typename MODEL::ObsAuxControl       ObsAuxControl_;
  typedef typename MODEL::ObsSpace            ObsSpace_;
  typedef typename MODEL::ObsVector           ObsVector_;
  typedef typename MODEL::Variables           Variables_;

 public:
  ObsOperatorBase() {}
  virtual ~ObsOperatorBase() {}

/// Obs Operator
  virtual void obsEquiv(const ModelAtLocations_ &, ObsVector_ &, const ObsAuxControl_ &) const =0;

/// Other
  virtual boost::shared_ptr<const Variables_> variables() const =0;  // Required from Model

 private:
  virtual void print(std::ostream &) const =0;
};

// =============================================================================

/// Obs Operator Factory
template <typename MODEL>
class ObsOperatorFactory {
  typedef typename MODEL::ObsSpace ObsSpace_;
 public:
  static ObsOperatorBase<MODEL> * create(const ObsSpace_ &, const eckit::Configuration &);
  virtual ~ObsOperatorFactory() { getMakers().clear(); }
 protected:
  explicit ObsOperatorFactory(const std::string &);
 private:
  virtual ObsOperatorBase<MODEL> * make(const ObsSpace_ &, const eckit::Configuration &) =0;
  static std::map < std::string, ObsOperatorFactory<MODEL> * > & getMakers() {
    static std::map < std::string, ObsOperatorFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class ObsOperatorMaker : public ObsOperatorFactory<MODEL> {
  typedef typename MODEL::ObsSpace const ObsSpace_;
  virtual ObsOperatorBase<MODEL> * make(const ObsSpace_ & odb, const eckit::Configuration & conf)
    { return new T(odb, conf); }
 public:
  explicit ObsOperatorMaker(const std::string & name) : ObsOperatorFactory<MODEL>(name) {}
};

// =============================================================================

template <typename MODEL>
ObsOperatorFactory<MODEL>::ObsOperatorFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    Log::error() << name << " already registered in observation operator factory." << std::endl;
    ABORT("Element already registered in ObsOperatorFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsOperatorBase<MODEL>* ObsOperatorFactory<MODEL>::create(const ObsSpace_ & odb,
                                                          const eckit::Configuration & conf) {
  Log::trace() << "ObsOperatorBase<MODEL>::create starting" << std::endl;
  Log::debug() << "ObsOperatorBase<MODEL>::create conf" << conf << std::endl;
  const std::string id = conf.getString("ObsType");
  typename std::map<std::string, ObsOperatorFactory<MODEL>*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    Log::error() << id << " does not exist in observation operator factory." << std::endl;
    ABORT("Element does not exist in ObsOperatorFactory.");
  }
  ObsOperatorBase<MODEL> * ptr = jloc->second->make(odb, conf);
  Log::trace() << "ObsOperatorBase<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSOPERATORBASE_H_
