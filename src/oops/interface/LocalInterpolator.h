/*
 * (C) Copyright 2022-2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/atlas/Interpolator.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"


namespace oops {
// Note: unlike many other "oops/interface/x" classes, the LocalInterpolator
// is NOT in the oops::interface namespace. That's because there is no
// derived "oops/base/x" class that wraps the interface for use by oops
// applications.

/// \brief Encapsulates local (ie on current MPI task) interpolators
template <typename MODEL>
class LocalInterpolator : private util::Printable,
                          private util::ObjectCounter<LocalInterpolator<MODEL>> {
  typedef typename MODEL::LocalInterpolator   LocalInterpolator_;
  typedef oops::Geometry<MODEL>            Geometry_;
  typedef oops::Increment<MODEL>           Increment_;
  typedef oops::State<MODEL>               State_;

  static constexpr bool IsGenericInterpolator =
    std::is_base_of_v<atlasbase::Interpolator, LocalInterpolator_>;

 public:
  static const std::string classname() {return "oops::LocalInterpolator";}

  LocalInterpolator(const eckit::Configuration &, const Geometry_ &,
                    const std::vector<double> &, const std::vector<double> &);
  ~LocalInterpolator();

  static void preprocess(State_ &);
  static void preprocess(Increment_ &);
  static void preprocessAD(Increment_ &);

  void apply(const Variables &, const State_ &,
             const std::vector<bool> &, std::vector<double> &) const;
  void apply(const Variables &, const Increment_ &,
             const std::vector<bool> &, std::vector<double> &) const;
  void applyAD(const Variables &, Increment_ &,
               const std::vector<bool> &, const std::vector<double> &) const;

 private:
  void print(std::ostream &) const override;

  std::unique_ptr<LocalInterpolator_> interpolator_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
LocalInterpolator<MODEL>::LocalInterpolator(const eckit::Configuration & conf,
                                            const Geometry_ & geometry,
                                            const std::vector<double> & lats,
                                            const std::vector<double> & lons)
  : interpolator_()
{
  Log::trace() << "LocalInterpolator<MODEL>::LocalInterpolator starting" << std::endl;
  util::Timer timer(classname(), "LocalInterpolator");
  if constexpr (IsGenericInterpolator) {
    interpolator_.reset(new LocalInterpolator_(conf, geometry.generic(), lats, lons));
  } else {
    interpolator_.reset(new LocalInterpolator_(conf, geometry.geometry(), lats, lons));
  }
  Log::trace() << "LocalInterpolator<MODEL>::LocalInterpolator done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
LocalInterpolator<MODEL>::~LocalInterpolator() {
  Log::trace() << "LocalInterpolator<MODEL>::~LocalInterpolator starting" << std::endl;
  util::Timer timer(classname(), "~LocalInterpolator");
  interpolator_.reset();
  Log::trace() << "LocalInterpolator<MODEL>::~LocalInterpolator done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalInterpolator<MODEL>::preprocess(State_ & xx) {
  Log::trace() << "LocalInterpolator<MODEL>::preprocess starting" << std::endl;
  util::Timer timer(classname(), "preprocess");
  if constexpr (IsGenericInterpolator) {
    LocalInterpolator_::preprocess(xx.fieldSet().fieldSet());
  } else {
    LocalInterpolator_::preprocess(xx.state());
  }
  Log::trace() << "LocalInterpolator<MODEL>::preprocess done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalInterpolator<MODEL>::preprocess(Increment_ & dx) {
  Log::trace() << "LocalInterpolator<MODEL>::preprocess starting" << std::endl;
  util::Timer timer(classname(), "preprocess");
  if constexpr (IsGenericInterpolator) {
    LocalInterpolator_::preprocess(dx.fieldSet().fieldSet());
  } else {
    LocalInterpolator_::preprocess(dx.increment());
  }
  Log::trace() << "LocalInterpolator<MODEL>::preprocess done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalInterpolator<MODEL>::preprocessAD(Increment_ & dx) {
  Log::trace() << "LocalInterpolator<MODEL>::preprocessAD starting" << std::endl;
  util::Timer timer(classname(), "preprocessAD");
  if constexpr (IsGenericInterpolator) {
    LocalInterpolator_::preprocessAD(dx.fieldSet().fieldSet());
    // Ensure data is propagated from the fieldset representation to the model-
    // specific representation at the end of the (adjoint) interpolation:
    dx.synchronizeFields();
  } else {
    LocalInterpolator_::preprocessAD(dx.increment());
  }
  Log::trace() << "LocalInterpolator<MODEL>::preprocessAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalInterpolator<MODEL>::apply(const Variables & vars, const State_ & xx,
                                     const std::vector<bool> & mask,
                                     std::vector<double> & buffer) const {
  Log::trace() << "LocalInterpolator<MODEL>::apply starting" << std::endl;
  util::Timer timer(classname(), "apply");
  if constexpr (IsGenericInterpolator) {
    interpolator_->apply(vars, xx.fieldSet().fieldSet(), mask, buffer);
  } else {
    interpolator_->apply(vars, xx.state(), mask, buffer);
  }
  Log::trace() << "LocalInterpolator<MODEL>::apply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalInterpolator<MODEL>::apply(const Variables & vars, const Increment_ & dx,
                                     const std::vector<bool> & mask,
                                     std::vector<double> & buffer) const {
  Log::trace() << "LocalInterpolator<MODEL>::applyTL starting" << std::endl;
  util::Timer timer(classname(), "applyTL");
  if constexpr (IsGenericInterpolator) {
    interpolator_->apply(vars, dx.fieldSet().fieldSet(), mask, buffer);
  } else {
    interpolator_->apply(vars, dx.increment(), mask, buffer);
  }
  Log::trace() << "LocalInterpolator<MODEL>::applyTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalInterpolator<MODEL>::applyAD(const Variables & vars, Increment_ & dx,
                                       const std::vector<bool> & mask,
                                       const std::vector<double> & buffer) const {
  Log::trace() << "LocalInterpolator<MODEL>::applyAD starting" << std::endl;
  util::Timer timer(classname(), "applyAD");
  if constexpr (IsGenericInterpolator) {
    interpolator_->applyAD(vars, dx.fieldSet().fieldSet(), mask, buffer);
  } else {
    interpolator_->applyAD(vars, dx.increment(), mask, buffer);
  }
  Log::trace() << "LocalInterpolator<MODEL>::applyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalInterpolator<MODEL>::print(std::ostream & os) const {
  Log::trace() << "LocalInterpolator<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *interpolator_;
  Log::trace() << "LocalInterpolator<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops
