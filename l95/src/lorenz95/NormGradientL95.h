/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef LORENZ95_NORMGRADIENTL95_H_
#define LORENZ95_NORMGRADIENTL95_H_

#include <memory>
#include <ostream>
#include <string>

#include "lorenz95/IncrementL95.h"
#include "lorenz95/Resolution.h"
#include "lorenz95/StateL95.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace lorenz95 {
  class Resolution;
  class IncrementL95;
  class StateL95;

// -----------------------------------------------------------------------------

class NormGradientL95 : public util::Printable,
                       private util::ObjectCounter<NormGradientL95> {
 public:
  static const std::string classname() {return "lorenz95::NormGradientL95";}

  NormGradientL95(const Resolution &, const StateL95 &, const eckit::Configuration &) {}
  ~NormGradientL95() {}

  void apply(IncrementL95 &) const {}

// Private
 private:
  void print(std::ostream & os) const {os << " NormGradientL95: print not implemented yet.";}
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_NORMGRADIENTL95_H_
