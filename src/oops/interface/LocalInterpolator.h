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
  interpolator_->apply(vars, xx.state(), mask, vect);
  Log::trace() << "LocalInterpolator<MODEL>::apply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalInterpolator<MODEL>::apply(const Variables & vars, const Increment_ & dx,
                                     const std::vector<bool> & mask,
                                     std::vector<double> & vect) const {
  Log::trace() << "LocalInterpolator<MODEL>::applyTL starting" << std::endl;
  util::Timer timer(classname(), "applyTL");
  interpolator_->apply(vars, dx.increment(), mask, vect);
  Log::trace() << "LocalInterpolator<MODEL>::applyTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalInterpolator<MODEL>::applyAD(const Variables & vars, Increment_ & dx,
                                       const std::vector<bool> & mask,
                                       const std::vector<double> & vect) const {
  Log::trace() << "LocalInterpolator<MODEL>::applyAD starting" << std::endl;
  util::Timer timer(classname(), "applyAD");
  interpolator_->applyAD(vars, dx.increment(), mask, vect);
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
