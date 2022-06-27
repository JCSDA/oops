/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"
#include "oops/util/TypeTraits.h"

namespace detail {

// ApplyHelper selects whether to call the model interpolator's apply and applyAD methods using an
// atlas::FieldSet interface or a model State/Increment interface.
//
// The fallback case uses the model State/Increment interface.
template <typename MODEL, typename = cpp17::void_t<>>
struct ApplyHelper {
  static void apply(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, const oops::State<MODEL> & state,
             const std::vector<bool> & mask, std::vector<double> & buffer) {
    interp.apply(vars, state.state(), mask, buffer);
  }
  static void apply(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, const oops::Increment<MODEL> & increment,
             const std::vector<bool> & mask, std::vector<double> & buffer) {
    interp.apply(vars, increment.increment(), mask, buffer);
  }
  static void applyAD(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, oops::Increment<MODEL> & increment,
             const std::vector<bool> & mask, const std::vector<double> & buffer) {
    interp.applyAD(vars, increment.increment(), mask, buffer);
  }
};

// The specialization calls the interpolator's FieldSet interface.
// Note: Here we *assume* that if the interpolator has apply(FieldSet interface), then it also
//       has applyAD(FieldSet interface). Code may fail to compile if this assumption is false.
template <typename MODEL>
struct ApplyHelper<MODEL, cpp17::void_t<decltype(std::declval<typename MODEL::LocalInterpolator>().
                                                   apply(std::declval<oops::Variables>(),
                                                         std::declval<atlas::FieldSet>(),
                                                         std::declval<std::vector<bool>>(),
                                                         std::declval<std::vector<double>&>()))>> {
  static void apply(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, const oops::State<MODEL> & state,
             const std::vector<bool> & mask, std::vector<double> & buffer) {
    interp.apply(vars, state.fieldSet(), mask, buffer);
  }
  static void apply(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, const oops::Increment<MODEL> & increment,
             const std::vector<bool> & mask, std::vector<double> & buffer) {
    interp.apply(vars, increment.fieldSet(), mask, buffer);
  }
  static void applyAD(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, oops::Increment<MODEL> & increment,
             const std::vector<bool> & mask, const std::vector<double> & buffer) {
    interp.applyAD(vars, increment.fieldSet(), mask, buffer);
  }
};
}  // namespace detail

namespace oops {

/// \brief Encapsulates local (ie on current MPI task) interpolators
// -----------------------------------------------------------------------------

template <typename MODEL>
class LocalInterpolator : public util::Printable,
                          private util::ObjectCounter<LocalInterpolator<MODEL> > {
  typedef typename MODEL::LocalInterpolator   LocalInterpolator_;
  typedef oops::Geometry<MODEL>            Geometry_;
  typedef oops::Increment<MODEL>           Increment_;
  typedef oops::State<MODEL>               State_;

 public:
  static const std::string classname() {return "oops::LocalInterpolator";}

  LocalInterpolator(const eckit::Configuration &, const Geometry_ &,
                    const std::vector<double> &, const std::vector<double> &);
  ~LocalInterpolator();

  void apply(const Variables &, const State_ &,
             const std::vector<bool> &, std::vector<double> &) const;
  void apply(const Variables &, const Increment_ &,
             const std::vector<bool> &, std::vector<double> &) const;
  void applyAD(const Variables &, Increment_ &,
               const std::vector<bool> &, const std::vector<double> &) const;

 private:
  std::unique_ptr<LocalInterpolator_> interpolator_;
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
LocalInterpolator<MODEL>::LocalInterpolator(const eckit::Configuration & conf,
                                            const Geometry_ & resol,
                                            const std::vector<double> & lats,
                                            const std::vector<double> & lons)
  : interpolator_()
{
  Log::trace() << "LocalInterpolator<MODEL>::LocalInterpolator starting" << std::endl;
  util::Timer timer(classname(), "LocalInterpolator");
  interpolator_.reset(new LocalInterpolator_(conf, resol.geometry(), lats, lons));
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
void LocalInterpolator<MODEL>::apply(const Variables & vars, const State_ & xx,
                                     const std::vector<bool> & mask,
                                     std::vector<double> & vect) const {
  Log::trace() << "LocalInterpolator<MODEL>::apply starting" << std::endl;
  util::Timer timer(classname(), "apply");
  detail::ApplyHelper<MODEL>::apply(*interpolator_, vars, xx, mask, vect);
  Log::trace() << "LocalInterpolator<MODEL>::apply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalInterpolator<MODEL>::apply(const Variables & vars, const Increment_ & dx,
                                     const std::vector<bool> & mask,
                                     std::vector<double> & vect) const {
  Log::trace() << "LocalInterpolator<MODEL>::applyTL starting" << std::endl;
  util::Timer timer(classname(), "applyTL");
  detail::ApplyHelper<MODEL>::apply(*interpolator_, vars, dx, mask, vect);
  Log::trace() << "LocalInterpolator<MODEL>::applyTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalInterpolator<MODEL>::applyAD(const Variables & vars, Increment_ & dx,
                                       const std::vector<bool> & mask,
                                       const std::vector<double> & vect) const {
  Log::trace() << "LocalInterpolator<MODEL>::applyAD starting" << std::endl;
  util::Timer timer(classname(), "applyAD");
  detail::ApplyHelper<MODEL>::applyAD(*interpolator_, vars, dx, mask, vect);
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
