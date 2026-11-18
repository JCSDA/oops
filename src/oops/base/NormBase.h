/*
 * (C) Crown Copyright 2026, Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_NORMBASE_H_
#define OOPS_BASE_NORMBASE_H_

#include <memory>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/Factory.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

// =============================================================================

/// \brief Virtual class representing a Norm object.
template <typename MODEL>
class NormBase : public util::Printable,
                 private util::ObjectCounter<NormBase<MODEL> > {
  typedef State<MODEL>            State_;
  typedef Increment<MODEL>        Increment_;

 public:
  virtual ~NormBase() = default;

  static const std::string classname() {return "oops::NormBase";}

// Apply norm to increment
  virtual void computeGradient(Increment_ &) const = 0;
  virtual void computeGradientInverse(Increment_ &) const = 0;
  virtual const double computeNorm(const Increment_ &) const = 0;
};

// =============================================================================

// Define a factory.
template<typename MODEL>
using NormFactory = util::Factory<NormBase<MODEL>,
                                  // Constructor args.
                                  const Variables &,
                                  const eckit::Configuration &,
                                  const State<MODEL> &>;

// =============================================================================

}  // namespace oops

#endif  // OOPS_BASE_NORMBASE_H_
