/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_INTERPOLATORBASE_H_
#define OOPS_BASE_INTERPOLATORBASE_H_

#include <map>
#include <memory>
#include <string>
#include <boost/noncopyable.hpp>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

/*! \brief Base class for Generic interpolation
 *
 * \details the **oops::InterpolatorBase** class provides a generic framework
 * for interpolation.
 *
 */

class InterpolatorBase : public util::Printable,
                         private boost::noncopyable {
 public:
  virtual ~InterpolatorBase() {}

  virtual void apply(atlas::FieldSet const &, atlas::FieldSet &) = 0;
  virtual void apply_ad(atlas::FieldSet const &, atlas::FieldSet &) = 0;

 private:
  virtual void print(std::ostream &) const = 0;
};

// -----------------------------------------------------------------------------
/// InterpolatorFactory: Factory for creating Interpolator objects

class InterpolatorFactory {
 public:
  static InterpolatorBase * create(const eckit::Configuration &,
                              const atlas::FunctionSpace &,
                              const atlas::FunctionSpace &,
                              const atlas::field::FieldSetImpl * = nullptr);
  virtual ~InterpolatorFactory() = default;
 protected:
  explicit InterpolatorFactory(const std::string &);
 private:
  virtual InterpolatorBase * make(const eckit::Configuration &,
                                  const atlas::FunctionSpace &,
                                  const atlas::FunctionSpace &,
                                  const atlas::field::FieldSetImpl *) = 0;
  static std::map < std::string, InterpolatorFactory * > & getMakers() {
    static std::map < std::string, InterpolatorFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------
/// Maker class for Interpolator Factory
///
/// Used to instantiate concrete Interpolator objects of
/// derived type T

template<class T>
class InterpolatorMaker : public InterpolatorFactory {
  virtual InterpolatorBase * make(const eckit::Configuration & conf,
                                  const atlas::FunctionSpace & fs1,
                                  const atlas::FunctionSpace & fs2,
                                  const atlas::field::FieldSetImpl * masks)
    { return new T(conf, fs1, fs2, masks); }
 public:
  explicit InterpolatorMaker(const std::string & name)
    : InterpolatorFactory(name) {}
};

// -----------------------------------------------------------------------------
/// Constructor for Interpolator Factory

InterpolatorFactory::InterpolatorFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    //
    // This was needed to get the bump interpolator to work with the gnu compilers
    // If the interpolator is already registered, do not abort.  Instead, just
    // write this message and return.
    //
    // Log::error() << name << " already registered in the interpolator factory."  << std::endl;
    // ABORT("Element already registered in InterpolatorFactory.");
    Log::info() << name << " already registered in the interpolator factory."  << std::endl;
  } else {
    getMakers()[name] = this;
  }
}

// -----------------------------------------------------------------------------
/// Create method for Interpolator Factory
///
/// This is what the user/developer will use to create Interpolator objects.
/// The default is to use atlas interpolation.

InterpolatorBase * InterpolatorFactory::create(
                                  const eckit::Configuration & conf,
                                  const atlas::FunctionSpace & fs1,
                                  const atlas::FunctionSpace & fs2,
                                  const atlas::field::FieldSetImpl * masks)
{
  Log::trace() << "InterpolatorBase::create starting" << std::endl;
  std::string id = conf.getString("interpolator", "atlas");
  typename std::map<std::string, InterpolatorFactory*>::iterator
    jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    Log::error() << id << " does not exist in the interpolator factory makers list." << std::endl;
    ABORT("Element does not exist in InterpolatorFactory.");
  }
  InterpolatorBase * ptr = jerr->second->make(conf, fs1, fs2, masks);
  Log::trace() << "InterpolatorBase::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_INTERPOLATORBASE_H_
