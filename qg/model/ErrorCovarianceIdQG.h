/*
 * (C) Copyright 2025- UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "oops/qg/GeometryQG.h"
#include "oops/qg/QgFortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace qg {
  class IncrementQG;
  class StateQG;

// -----------------------------------------------------------------------------
/// Identity background error covariance matrix for QG model.

class ErrorCovarianceIdQG : public util::Printable,
                            private boost::noncopyable,
                            private util::ObjectCounter<ErrorCovarianceIdQG> {
 public:
  static const std::string classname() {return "qg::ErrorCovarianceIdQG";}

  ErrorCovarianceIdQG(const GeometryQG &, const oops::Variables &,
                      const eckit::Configuration &, const StateQG &,
                      const StateQG &);
  ~ErrorCovarianceIdQG() = default;

  void multiply(const IncrementQG &, IncrementQG &) const;
  void inverseMultiply(const IncrementQG &, IncrementQG &) const;
  void randomize(IncrementQG &) const;

 private:
  void print(std::ostream &) const;
};
// -----------------------------------------------------------------------------

}  // namespace qg
