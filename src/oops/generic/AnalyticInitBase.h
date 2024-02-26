/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_ANALYTICINITBASE_H_
#define OOPS_GENERIC_ANALYTICINITBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "oops/interface/GeoVaLs.h"
#include "oops/interface/SampledLocations.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/ObjectCounter.h"

namespace oops {

/// \brief Initializes GeoVaLs with analytic formula
template <typename OBS>
class AnalyticInitBase : private util::ObjectCounter<AnalyticInitBase<OBS> > {
  typedef GeoVaLs<OBS>                GeoVaLs_;
  typedef SampledLocations<OBS>       SampledLocations_;

 public:
  static const std::string classname() {return "oops::AnalyticInitBase";}

  AnalyticInitBase() = default;
  virtual ~AnalyticInitBase() = default;

/*! \brief Fill GeoVaLs with values computed by analytic function.
 *
 * \details **AnalyticInit::fillGeoVaLs** was introduced in May, 2018 (initially
 * as a GeoVaLs constructor) for use with the interpolation test. The interpolation test
 * requires an initialization of a GeoVaLs object based on the same analytic
 * formulae used for the State initialization (see test::TestStateInterpolation()
 * for further information).  This in turn requires information about the
 * vertical profile in addition to the latitude and longitude positional
 * information in the SampledLocations object.  Currently, this information
 * about the vertical profile is obtained from an existing GeoVaLs object
 * (passed as *gvals*) that represents the output of the State::interpolate()
 * method.
 *
 * \date May, 2018: created as a constructor (M. Miesch, JCSDA)
 * \date June, 2018: moved to a method (M. Miesch, JCSDA)
 *
 * \sa test::TestStateInterpolation()
 */
  virtual void fillGeoVaLs(const SampledLocations_ &, GeoVaLs_ &) const = 0;
};

// -----------------------------------------------------------------------------

/// A factory creating instances of concrete subclasses of AnalyticInitBase.
template <typename OBS>
class AnalyticInitFactory {
  typedef AnalyticInitBase<OBS> AnalyticInitBase_;

 public:
  static std::unique_ptr<AnalyticInitBase_> create(const eckit::Configuration &);

  /// \brief Return the names of all methods that can be created by one of the registered makers.
  static std::vector<std::string> getMakerNames() {
    return keys(getMakers());
  }

  virtual ~AnalyticInitFactory() = default;

 protected:
  explicit AnalyticInitFactory(const std::string &);

 private:
  virtual std::unique_ptr<AnalyticInitBase_> make(const eckit::Configuration &) = 0;

  static std::map < std::string, AnalyticInitFactory<OBS> * > & getMakers() {
    static std::map < std::string, AnalyticInitFactory<OBS> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

/// \brief A subclass of AnalyticInitFactory able to create instances of T (a concrete subclass of
/// AnalyticInitBase<OBS>).
template<typename OBS, typename T>
class AnalyticInitMaker : public AnalyticInitFactory<OBS> {
  typedef AnalyticInitBase<OBS>       AnalyticInitBase_;

  std::unique_ptr<AnalyticInitBase_> make(const eckit::Configuration & conf) override {
    return std::make_unique<T>(conf);
  }

 public:
  explicit AnalyticInitMaker(const std::string & name) : AnalyticInitFactory<OBS>(name) {}
};

// -----------------------------------------------------------------------------

template <typename OBS>
AnalyticInitFactory<OBS>::AnalyticInitFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in analytic init factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename OBS>
std::unique_ptr<AnalyticInitBase<OBS>>
AnalyticInitFactory<OBS>::create(const eckit::Configuration & config) {
  Log::trace() << "AnalyticInitFactory<OBS>::create starting" << std::endl;
  const std::string &id = config.getString("method");
  typename std::map<std::string, AnalyticInitFactory<OBS>*>::iterator
    jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    throw std::runtime_error(id + " does not exist in analytic init factory.");
  }
  std::unique_ptr<AnalyticInitBase<OBS>> ptr(jerr->second->make(config));
  Log::trace() << "AnalyticInitFactory<OBS>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_ANALYTICINITBASE_H_
