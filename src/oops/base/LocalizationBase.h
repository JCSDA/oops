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

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/mpi/mpi.h"
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

 public:
  LocalizationBase() = default;
  virtual ~LocalizationBase() = default;

  /// 4D localization with the same localization for time blocks
  virtual void randomize(Increment_ &) const;
  virtual void multiply(Increment_ &) const;

 protected:
  virtual void doRandomize(Increment_ &) const = 0;
  virtual void doMultiply(Increment_ &) const = 0;

 private:
  virtual void print(std::ostream &) const = 0;
};

// -----------------------------------------------------------------------------
/// Randomize
template <typename MODEL>
void LocalizationBase<MODEL>::randomize(Increment_ & dx) const {
  Log::trace() << "LocalizationBase<MODEL>::randomize starting" << std::endl;
  const eckit::mpi::Comm & comm = dx.timeComm();
  static int tag = 23456;
  size_t nslots = comm.size();
  int mytime = comm.rank();

  if (mytime > 0) {
    util::DateTime dt = dx.validTime();   // Save original time value
    dx.zero();
    oops::mpi::receive(comm, dx, 0, tag);
    dx.updateTime(dt - dx.validTime());  // Set time back to original value
  } else {
    // Apply 3D localization
    this->doRandomize(dx);

    // Copy result to all timeslots
    for (size_t jj = 1; jj < nslots; ++jj) {
      oops::mpi::send(comm, dx, jj, tag);
    }
  }
  ++tag;

  Log::trace() << "LocalizationBase<MODEL>::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------
/// Multiply
template <typename MODEL>
void LocalizationBase<MODEL>::multiply(Increment_ & dx) const {
  Log::trace() << "LocalizationBase<MODEL>::multiply starting" << std::endl;
  const eckit::mpi::Comm & comm = dx.timeComm();
  static int tag = 23456;
  size_t nslots = comm.size();
  int mytime = comm.rank();

  if (mytime > 0) {
    util::DateTime dt = dx.validTime();   // Save original time value
    oops::mpi::send(comm, dx, 0, tag);
    dx.zero();
    oops::mpi::receive(comm, dx, 0, tag);
    dx.updateTime(dt - dx.validTime());  // Set time back to original value
  } else {
    // Sum over timeslots
    for (size_t jj = 1; jj < nslots; ++jj) {
      Increment_ dxtmp(dx);
      oops::mpi::receive(comm, dxtmp, jj, tag);
      dx.axpy(1.0, dxtmp, false);
    }

    // Apply 3D localization
    this->doMultiply(dx);

    // Copy result to all timeslots
    for (size_t jj = 1; jj < nslots; ++jj) {
      oops::mpi::send(comm, dx, jj, tag);
    }
  }
  ++tag;

  Log::trace() << "LocalizationBase<MODEL>::multiply done" << std::endl;
}

// =============================================================================

/// Localization Factory
template <typename MODEL>
class LocalizationFactory {
  typedef Geometry<MODEL>                             Geometry_;
 public:
  static std::unique_ptr<LocalizationBase<MODEL>> create(const Geometry_ &,
                                                         const util::DateTime &,
                                                         const eckit::Configuration &);
  virtual ~LocalizationFactory() = default;
 protected:
  explicit LocalizationFactory(const std::string &);
 private:
  virtual LocalizationBase<MODEL> * make(const Geometry_ &,
                                         const util::DateTime &,
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
                                         const util::DateTime & time,
                                         const eckit::Configuration & conf)
    { return new T(geometry, time, conf); }
 public:
  explicit LocalizationMaker(const std::string & name) : LocalizationFactory<MODEL>(name) {}
};

// =============================================================================

template <typename MODEL>
LocalizationFactory<MODEL>::LocalizationFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in localization factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
std::unique_ptr<LocalizationBase<MODEL>>
LocalizationFactory<MODEL>::create(const Geometry_ & geometry,
                                   const util::DateTime & time,
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
    throw std::runtime_error(id + " does not exist in localization factory.");
  }
  std::unique_ptr<LocalizationBase<MODEL>> ptr(jloc->second->make(geometry, time, conf));
  Log::trace() << "LocalizationBase<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_LOCALIZATIONBASE_H_
