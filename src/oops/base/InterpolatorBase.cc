/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/InterpolatorBase.h"

#include <string>


#include "eckit/config/Configuration.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace oops {
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
