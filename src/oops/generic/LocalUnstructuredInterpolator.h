/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/generic/InterpolatorUnstructured.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------

template<typename MODEL>
class LocalUnstructuredInterpolator : public util::Printable {
  typedef State<MODEL>     State_;
  typedef Geometry<MODEL>  Geometry_;
  typedef Increment<MODEL> Increment_;
 public:
  LocalUnstructuredInterpolator(const eckit::Configuration &, const Geometry_ &,
                                const std::vector<double> &);

  void apply(const Variables &, const State_ &, std::vector<double> &) const;
  void apply(const Variables &, const Increment_ &, std::vector<double> &) const;
  void applyAD(const Variables &, Increment_ &, const std::vector<double> &) const;

 private:
  void print(std::ostream &) const;

  std::unique_ptr<InterpolatorUnstructured> interp_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
LocalUnstructuredInterpolator<MODEL>::
  LocalUnstructuredInterpolator(const eckit::Configuration & config, const Geometry_ & grid,
                                const std::vector<double> & latlon_out) : interp_()
{
  Log::trace() << "LocalUnstructuredInterpolator::LocalUnstructuredInterpolator start" << std::endl;

  std::vector<double> lats_in;
  std::vector<double> lons_in;
  grid.latlon(lats_in, lons_in, true);

  interp_ = std::make_unique<InterpolatorUnstructured>(config, lats_in, lons_in, latlon_out);

  Log::trace() << "LocalUnstructuredInterpolator::LocalUnstructuredInterpolator done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalUnstructuredInterpolator<MODEL>::
  apply(const Variables & vars, const State_ & xx, std::vector<double> & locvals) const
{
  interp_->apply(vars, xx.fieldSet(), locvals);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalUnstructuredInterpolator<MODEL>::
  apply(const Variables & vars, const Increment_ & dx, std::vector<double> & locvals) const
{
  interp_->apply(vars, dx.fieldSet(), locvals);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalUnstructuredInterpolator<MODEL>::
  applyAD(const Variables & vars, Increment_ & dx, const std::vector<double> & locvals) const
{
  interp_->applyAD(vars, dx.fieldSet(), locvals);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalUnstructuredInterpolator<MODEL>::print(std::ostream & os) const
{
  os << "LocalUnstructuredInterpolator<" << MODEL::name() << ">";
}

// -----------------------------------------------------------------------------

}  // namespace oops
