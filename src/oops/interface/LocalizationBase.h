/*
* Copyright 2011 ECMWF
*
* This software was developed at ECMWF for evaluation
* and may be used for academic and research purposes only.
* The software is provided as is without any warranty.
*
* This software can be used, copied and modified but not
* redistributed or sold. This notice must be reproduced
* on each copy made.
*/

#ifndef OOPS_INTERFACE_LOCALIZATIONBASE_H_
#define OOPS_INTERFACE_LOCALIZATIONBASE_H_

#include <mpi.h>
#include <map>
#include <string>
#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "eckit/config/Configuration.h"
#include "oops/assimilation/Increment4D.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Base class for localizations

template<typename MODEL>
class LocalizationBase : public util::Printable,
                         private boost::noncopyable {
  typedef Increment<MODEL>           Increment_;
  typedef Increment4D<MODEL>         Increment4D_;

 public:
  LocalizationBase() {}
  virtual ~LocalizationBase() {}

  void doMultiply(Increment_ &) const;
  void doMultiply(Increment4D_ &) const;

 private:
  virtual void multiply(typename MODEL::Increment &) const = 0;
  virtual void print(std::ostream &) const = 0;
};

// =============================================================================

/// LocalizationFactory Factory
template <typename MODEL>
class LocalizationFactory {
  typedef Geometry<MODEL> Geometry_;
 public:
  static LocalizationBase<MODEL> * create(const Geometry_ &, const eckit::Configuration &);
  virtual ~LocalizationFactory() { getMakers().clear(); }
 protected:
  explicit LocalizationFactory(const std::string &);
 private:
  virtual LocalizationBase<MODEL> * make(const Geometry_ &, const eckit::Configuration &) = 0;
  static std::map < std::string, LocalizationFactory<MODEL> * > & getMakers() {
    static std::map < std::string, LocalizationFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class LocalizationMaker : public LocalizationFactory<MODEL> {
  typedef Geometry<MODEL> Geometry_;
  virtual LocalizationBase<MODEL> * make(const Geometry_ & resol, const eckit::Configuration & conf)
    { return new T(resol.geometry(), conf); }
 public:
  explicit LocalizationMaker(const std::string & name) : LocalizationFactory<MODEL>(name) {}
};

// -----------------------------------------------------------------------------

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
LocalizationBase<MODEL>* LocalizationFactory<MODEL>::create(const Geometry_ & resol,
                                                            const eckit::Configuration & conf) {
  Log::trace() << "LocalizationBase<MODEL>::create starting" << std::endl;
  const std::string id = conf.getString("localization");
  typename std::map<std::string, LocalizationFactory<MODEL>*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    Log::trace() << id << " does not exist in localization factory." << std::endl;
    ABORT("Element does not exist in LocalizationFactory.");
  }
  LocalizationBase<MODEL> * ptr = jloc->second->make(resol, conf);
  Log::trace() << "LocalizationBase<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LocalizationBase<MODEL>::doMultiply(Increment_ & dx) const {
  Log::trace() << "LocalizationBase<MODEL>::multiply starting" << std::endl;
  this->multiply(dx.increment());
  Log::trace() << "LocalizationBase<MODEL>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LocalizationBase<MODEL>::doMultiply(Increment4D_ & dx) const {
  Log::trace() << "LocalizationBase<MODEL>::multiply starting" << std::endl;
  // Sum over timeslots
  Increment_ dxtmp(dx[dx.first()]);
  for (int isub = dx.first()+1; isub <= dx.last(); ++isub) {
     dxtmp.axpy(1.0, dx[isub], false);
  }

  // Apply 3D localization
  this->multiply(dxtmp.increment());

  // Copy result to all timeslots
  for (int isub = dx.first(); isub <= dx.last(); ++isub) {
     dx[isub].zero();
     dx[isub].axpy(1.0, dxtmp, false);
  }
  Log::trace() << "LocalizationBase<MODEL>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_LOCALIZATIONBASE_H_
