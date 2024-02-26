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

#include "atlas/field/FieldSet.h"
#include "eckit/config/Configuration.h"
#include "oops/base/Geometry.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/generic/LocalInterpolatorBase.h"
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
  // without mask
  static void apply(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, const oops::State<MODEL> & state,
             std::vector<double> & buffer) {
    interp.apply(vars, state.state(), buffer);
  }
  static void apply(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, const oops::Increment<MODEL> & increment,
             std::vector<double> & buffer) {
    interp.apply(vars, increment.increment(), buffer);
  }
  static void applyAD(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, oops::Increment<MODEL> & increment,
             const std::vector<double> & buffer) {
    interp.applyAD(vars, increment.increment(), buffer);
  }

  // with mask
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
// Note: Here we *assume* that...
// - if the interpolator has apply(FieldSet interface), then it has applyAD(FieldSet interface)
// - if the interpolator has apply(with mask), then it has apply(without mask)
// Code may fail to compile if either assumption is violated.
template <typename MODEL>
struct ApplyHelper<MODEL, cpp17::void_t<decltype(std::declval<typename MODEL::LocalInterpolator>().
                                                   apply(std::declval<oops::Variables>(),
                                                         std::declval<atlas::FieldSet>(),
                                                         std::declval<std::vector<bool>>(),
                                                         std::declval<std::vector<double>&>()))>> {
  // without mask
  static void apply(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, const oops::State<MODEL> & state,
             std::vector<double> & buffer) {
    interp.apply(vars, state.fieldSet().fieldSet(), buffer);
  }
  static void apply(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, const oops::Increment<MODEL> & increment,
             std::vector<double> & buffer) {
    interp.apply(vars, increment.fieldSet().fieldSet(), buffer);
  }
  static void applyAD(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, oops::Increment<MODEL> & increment,
             const std::vector<double> & buffer) {
    interp.applyAD(vars, increment.fieldSet().fieldSet(), buffer);
  }

  // with mask
  static void apply(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, const oops::State<MODEL> & state,
             const std::vector<bool> & mask, std::vector<double> & buffer) {
    interp.apply(vars, state.fieldSet().fieldSet(), mask, buffer);
  }
  static void apply(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, const oops::Increment<MODEL> & increment,
             const std::vector<bool> & mask, std::vector<double> & buffer) {
    interp.apply(vars, increment.fieldSet().fieldSet(), mask, buffer);
  }
  static void applyAD(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, oops::Increment<MODEL> & increment,
             const std::vector<bool> & mask, const std::vector<double> & buffer) {
    interp.applyAD(vars, increment.fieldSet().fieldSet(), mask, buffer);
  }
};

// Constructor helper - determine whether interpolator has a GeomData based constructor
template <typename MODEL, typename = void>
struct ConstructHelper {
  static typename MODEL::LocalInterpolator* apply(const eckit::Configuration & conf,
                                                  const oops::Geometry<MODEL> & geom,
                                                  const std::vector<double> & lats,
                                                  const std::vector<double> & lons) {
    return new typename MODEL::LocalInterpolator(conf, geom.geometry(), lats, lons);
  }
};

template <typename MODEL>
struct ConstructHelper<MODEL, cpp17::void_t<typename std::enable_if<std::is_constructible<
                              typename MODEL::LocalInterpolator,
                              const eckit::Configuration &,
                              const oops::GeometryData &,
                              const std::vector<double> &,
                              const std::vector<double> &>::value>::type>>{
  static typename MODEL::LocalInterpolator* apply(const eckit::Configuration & conf,
                                                  const oops::Geometry<MODEL> & geom,
                                                  const std::vector<double> & lats,
                                                  const std::vector<double> & lons) {
    return new typename MODEL::LocalInterpolator(conf, geom.generic(), lats, lons);
  }
};

// ApplyAtlasHelper tries to call the model interpolator's apply(AD) methods.
// The fallback case errors.
template <typename MODEL, typename = cpp17::void_t<>>
struct ApplyAtlasHelper {
  // without target-point mask
  static void apply(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, const atlas::FieldSet & fset,
             std::vector<double> & buffer) {
    throw eckit::Exception("The LocalInterpolator for model " + MODEL::name() + " has no method "
                           "apply (taking an atlas::FieldSet), but an oops::LocalInterpolator "
                           "tried to call this non-existent method.");
  }
  static void applyAD(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, atlas::FieldSet & fset,
             const std::vector<double> & buffer) {
    throw eckit::Exception("The LocalInterpolator for model " + MODEL::name() + " has no method "
                           "applyAD (taking an atlas::FieldSet), but an oops::LocalInterpolator "
                           "tried to call this non-existent method.");
  }
  // with target-point mask
  static void apply(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, const atlas::FieldSet & fset,
             const std::vector<bool> & mask, std::vector<double> & buffer) {
    throw eckit::Exception("The LocalInterpolator for model " + MODEL::name() + " has no method "
                           "apply (taking an atlas::FieldSet), but an oops::LocalInterpolator "
                           "tried to call this non-existent method.");
  }
  static void applyAD(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, atlas::FieldSet & fset,
             const std::vector<bool> & mask, const std::vector<double> & buffer) {
    throw eckit::Exception("The LocalInterpolator for model " + MODEL::name() + " has no method "
                           "applyAD (taking an atlas::FieldSet), but an oops::LocalInterpolator "
                           "tried to call this non-existent method.");
  }
};

// The specialization calls the model-specific interpolator's atlas::FieldSet interface.
// Note: Here we simplify the template metaprogramming by assuming that if the interpolator has a
//       method apply(FieldSet), then it also has applyAD(FieldSet), and also assuming that
//       interfaces without and with target-point masks exist. The code may fail to compile if these
//       assumptions are violated.
template <typename MODEL>
struct ApplyAtlasHelper<MODEL,
    cpp17::void_t<decltype(std::declval<typename MODEL::LocalInterpolator>().apply(
        std::declval<oops::Variables>(),
        std::declval<atlas::FieldSet>(),
        std::declval<std::vector<double>&>()))>> {
  // without target-point mask
  static void apply(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, const atlas::FieldSet & fset,
             std::vector<double> & buffer) {
    interp.apply(vars, fset, buffer);
  }
  static void applyAD(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, atlas::FieldSet & fset,
             const std::vector<double> & buffer) {
    interp.applyAD(vars, fset, buffer);
  }
  // with target-point mask
  static void apply(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, const atlas::FieldSet & fset,
             const std::vector<bool> & mask, std::vector<double> & buffer) {
    interp.apply(vars, fset, mask, buffer);
  }
  static void applyAD(const typename MODEL::LocalInterpolator & interp,
             const oops::Variables & vars, atlas::FieldSet & fset,
             const std::vector<bool> & mask, const std::vector<double> & buffer) {
    interp.applyAD(vars, fset, mask, buffer);
  }
};

}  // namespace detail

namespace oops {

/// \brief Encapsulates local (ie on current MPI task) interpolators
// -----------------------------------------------------------------------------

template <typename MODEL>
class LocalInterpolator : public LocalInterpolatorBase,
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
  void apply(const Variables &, const State_ &,
             std::vector<double> &) const;
  void apply(const Variables &, const Increment_ &,
             const std::vector<bool> &, std::vector<double> &) const;
  void apply(const Variables &, const Increment_ &,
             std::vector<double> &) const;
  void apply(const Variables &, const atlas::FieldSet &,
             const std::vector<bool> &, std::vector<double> &) const override;
  void apply(const Variables &, const atlas::FieldSet &,
             std::vector<double> &) const override;
  void applyAD(const Variables &, Increment_ &,
               const std::vector<bool> &, const std::vector<double> &) const;
  void applyAD(const Variables &, Increment_ &,
               const std::vector<double> &) const;
  void applyAD(const Variables &, atlas::FieldSet &,
               const std::vector<bool> &, const std::vector<double> &) const override;
  void applyAD(const Variables &, atlas::FieldSet &,
               const std::vector<double> &) const override;

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
  interpolator_.reset(::detail::ConstructHelper<MODEL>::apply(conf, resol, lats, lons));
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
  ::detail::ApplyHelper<MODEL>::apply(*interpolator_, vars, xx, mask, vect);
  Log::trace() << "LocalInterpolator<MODEL>::apply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalInterpolator<MODEL>::apply(const Variables & vars, const State_ & xx,
                                     std::vector<double> & vect) const {
  Log::trace() << "LocalInterpolator<MODEL>::apply starting" << std::endl;
  util::Timer timer(classname(), "apply");
  ::detail::ApplyHelper<MODEL>::apply(*interpolator_, vars, xx, vect);
  Log::trace() << "LocalInterpolator<MODEL>::apply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalInterpolator<MODEL>::apply(const Variables & vars, const Increment_ & dx,
                                     const std::vector<bool> & mask,
                                     std::vector<double> & vect) const {
  Log::trace() << "LocalInterpolator<MODEL>::applyTL starting" << std::endl;
  util::Timer timer(classname(), "applyTL");
  ::detail::ApplyHelper<MODEL>::apply(*interpolator_, vars, dx, mask, vect);
  Log::trace() << "LocalInterpolator<MODEL>::applyTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalInterpolator<MODEL>::apply(const Variables & vars, const Increment_ & dx,
                                     std::vector<double> & vect) const {
  Log::trace() << "LocalInterpolator<MODEL>::applyTL starting" << std::endl;
  util::Timer timer(classname(), "applyTL");
  ::detail::ApplyHelper<MODEL>::apply(*interpolator_, vars, dx, vect);
  Log::trace() << "LocalInterpolator<MODEL>::applyTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalInterpolator<MODEL>::apply(const Variables & vars, const atlas::FieldSet & fs,
                                     const std::vector<bool> & mask,
                                     std::vector<double> & vect) const {
  Log::trace() << "LocalInterpolator<MODEL>::apply(FieldSet) starting" << std::endl;
  util::Timer timer(classname(), "applyTL");
  ::detail::ApplyAtlasHelper<MODEL>::apply(*interpolator_, vars, fs, mask, vect);
  Log::trace() << "LocalInterpolator<MODEL>::apply(FieldSet) done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalInterpolator<MODEL>::apply(const Variables & vars, const atlas::FieldSet & fs,
                                     std::vector<double> & vect) const {
  Log::trace() << "LocalInterpolator<MODEL>::apply(FieldSet) starting" << std::endl;
  util::Timer timer(classname(), "applyTL");
  ::detail::ApplyAtlasHelper<MODEL>::apply(*interpolator_, vars, fs, vect);
  Log::trace() << "LocalInterpolator<MODEL>::apply(FieldSet) done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalInterpolator<MODEL>::applyAD(const Variables & vars, Increment_ & dx,
                                       const std::vector<bool> & mask,
                                       const std::vector<double> & vect) const {
  Log::trace() << "LocalInterpolator<MODEL>::applyAD starting" << std::endl;
  util::Timer timer(classname(), "applyAD");
  ::detail::ApplyHelper<MODEL>::applyAD(*interpolator_, vars, dx, mask, vect);
  Log::trace() << "LocalInterpolator<MODEL>::applyAD done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL>
void LocalInterpolator<MODEL>::applyAD(const Variables & vars, Increment_ & dx,
                                       const std::vector<double> & vect) const {
  Log::trace() << "LocalInterpolator<MODEL>::applyAD starting" << std::endl;
  util::Timer timer(classname(), "applyAD");
  ::detail::ApplyHelper<MODEL>::applyAD(*interpolator_, vars, dx, vect);
  Log::trace() << "LocalInterpolator<MODEL>::applyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalInterpolator<MODEL>::applyAD(const Variables & vars, atlas::FieldSet & fs,
                                       const std::vector<bool> & mask,
                                       const std::vector<double> & vect) const {
  Log::trace() << "LocalInterpolator<MODEL>::applyAD(FieldSet) starting" << std::endl;
  util::Timer timer(classname(), "applyAD");
  ::detail::ApplyAtlasHelper<MODEL>::applyAD(*interpolator_, vars, fs, mask, vect);
  Log::trace() << "LocalInterpolator<MODEL>::applyAD(FieldSet) done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalInterpolator<MODEL>::applyAD(const Variables & vars, atlas::FieldSet & fs,
                                       const std::vector<double> & vect) const {
  Log::trace() << "LocalInterpolator<MODEL>::applyAD(FieldSet) starting" << std::endl;
  util::Timer timer(classname(), "applyAD");
  ::detail::ApplyAtlasHelper<MODEL>::applyAD(*interpolator_, vars, fs, vect);
  Log::trace() << "LocalInterpolator<MODEL>::applyAD(FieldSet) done" << std::endl;
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
