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

#ifndef OOPS_GENERIC_LOCALIZATIONBASE_H_
#define OOPS_GENERIC_LOCALIZATIONBASE_H_

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
/// \brief Base class for generic implementations of model-space localization.
/// Use this class as a base class for generic implementations,
/// and interface::LocalizationBase as a base class for MODEL-specific implementations.
///
/// Note: generic implementations need to provide a constructor with the following signature:
///
///     LocalizationBase(const Geometry<MODEL> &,
///                      const oops::Variables &,
///                      const eckit::Configuration &);
template <typename MODEL>
class LocalizationBase : public util::Printable,
                         private boost::noncopyable {
  typedef Increment<MODEL>  Increment_;

 public:
  LocalizationBase() = default;
  virtual ~LocalizationBase() = default;

  /// Randomize \p dx and apply 3D localization
  virtual void randomize(Increment_ & dx) const = 0;
  /// Apply 3D localization to \p dx
  virtual void multiply(Increment_ & dx) const = 0;
};

// -----------------------------------------------------------------------------

/// Localization Factory
template <typename MODEL>
class LocalizationFactory {
  typedef Geometry<MODEL>                             Geometry_;
 public:
  /// \brief Create and return a new Localization.
  ///
  /// The Localization type is determined by the "localization method" entry of
  /// \p config.
  static std::unique_ptr<LocalizationBase<MODEL>> create(const Geometry_ &,
                                                         const oops::Variables &,
                                                         const eckit::Configuration & config);
  virtual ~LocalizationFactory() = default;
 protected:
  /// \brief Register a maker able to create Localizations of type \p name.
  explicit LocalizationFactory(const std::string & name);
 private:
  virtual std::unique_ptr<LocalizationBase<MODEL>> make(const Geometry_ &,
                                                        const oops::Variables &,
                                                        const eckit::Configuration &) = 0;
  static std::map < std::string, LocalizationFactory<MODEL> * > & getMakers() {
    static std::map < std::string, LocalizationFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------
/// \brief A subclass of LocalizationFactory able to create instances of T (a concrete
/// subclass of LocalizationBase<MODEL>). Passes Geometry<MODEL> to the constructor of T.
template<class MODEL, class T>
class LocalizationMaker : public LocalizationFactory<MODEL> {
  typedef Geometry<MODEL>                             Geometry_;
  std::unique_ptr<LocalizationBase<MODEL>> make(const Geometry_ & geometry,
                                                const oops::Variables & incVars,
                                                const eckit::Configuration & conf) override
    { return std::make_unique<T>(geometry, incVars, conf); }
 public:
  explicit LocalizationMaker(const std::string & name) : LocalizationFactory<MODEL>(name) {}
};

// -----------------------------------------------------------------------------

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
                                   const oops::Variables & incVars,
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
  std::unique_ptr<LocalizationBase<MODEL>> ptr(jloc->second->make(geometry, incVars, conf));
  Log::trace() << "LocalizationBase<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_LOCALIZATIONBASE_H_
