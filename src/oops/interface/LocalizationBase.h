/*
* Copyright 2011 ECMWF
* Copyright 2020-2021 UCAR
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

#include <memory>
#include <string>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/generic/LocalizationBase.h"

namespace eckit {
  class Configuration;
}

namespace oops {

namespace interface {

// -----------------------------------------------------------------------------
/// Base class for MODEL-specific implementations of the Loclaization interface.
/// interface::LocalizationBase overrides oops::LocalizationBase methods to pass
/// MODEL-specific implementation of Increment to the MODEL-specific
/// implementation of Localization.
///
/// Note: subclasses need to provide a constructor with the following signature:
///
///     LocalizationBase(const MODEL::Geometry &,
///                      const oops::Variables &,
///                      const eckit::Configuration &);
///
template<typename MODEL>
class LocalizationBase : public oops::LocalizationBase<MODEL> {
  typedef typename MODEL::Increment   Increment_;
 public:
  static const std::string classname() {return "oops::Localization";}

  LocalizationBase() = default;
  virtual ~LocalizationBase() = default;

  /// Overrides for oops::LocalizationBase classes, passing MODEL-specific classes to the
  /// MODEL-specific implementations of Localization
  void randomize(oops::Increment<MODEL> & dx) const final
       { this->randomize(dx.increment()); }
  void multiply(oops::Increment<MODEL> & dx) const final
       { this->multiply(dx.increment()); }

  /// Randomize \p dx and apply 3D localization
  virtual void randomize(Increment_ & dx) const = 0;
  /// Apply 3D localization to \p dx
  virtual void multiply(Increment_ & dx) const = 0;
};

// -----------------------------------------------------------------------------
/// \brief A subclass of LocalizationFactory able to create instances of T (a concrete
/// subclass of interface::LocalizationBase<MODEL>). Passes MODEL::Geometry to the constructor of T.
template<class MODEL, class T>
class LocalizationMaker : public oops::LocalizationFactory<MODEL> {
  typedef oops::Geometry<MODEL>  Geometry_;
  std::unique_ptr<oops::LocalizationBase<MODEL>> make(const Geometry_ & geometry,
                                                      const oops::Variables & vars,
                                                      const eckit::Configuration & conf) override
    { return std::make_unique<T>(geometry.geometry(), vars, conf); }
 public:
  explicit LocalizationMaker(const std::string & name) : oops::LocalizationFactory<MODEL>(name) {}
};


// -----------------------------------------------------------------------------

}  // namespace interface

}  // namespace oops

#endif  // OOPS_INTERFACE_LOCALIZATIONBASE_H_
