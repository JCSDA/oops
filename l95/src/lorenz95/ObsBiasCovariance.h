/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_OBSBIASCOVARIANCE_H_
#define LORENZ95_OBSBIASCOVARIANCE_H_

#include <memory>
#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "lorenz95/ObsBiasPreconditioner.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace lorenz95 {
  class ObsBias;
  class ObsBiasCorrection;
  class ObsTable;

// -----------------------------------------------------------------------------

class ObsBiasCovariance : public util::Printable,
                          private boost::noncopyable,
                          private util::ObjectCounter<ObsBiasCovariance> {
 public:
  static const std::string classname() {return "lorenz95::ObsBiasCovariance";}

/// Constructor, destructor
  ObsBiasCovariance(const ObsTable &, const eckit::Configuration &);
  ~ObsBiasCovariance() {}

/// Linear algebra operators
  void linearize(const ObsBias &, const eckit::Configuration &) {}
  void multiply(const ObsBiasCorrection &, ObsBiasCorrection &) const;
  void inverseMultiply(const ObsBiasCorrection &, ObsBiasCorrection &) const;
  void randomize(ObsBiasCorrection &) const;
  std::unique_ptr<ObsBiasPreconditioner> preconditioner() const;

  /// I/O and diagnostics
  void write(const eckit::Configuration &) const {}
  bool active() const {return active_;}

 private:
  void print(std::ostream &) const;
  double variance_;
  bool active_;
};

// -----------------------------------------------------------------------------
}  // namespace lorenz95

#endif  // LORENZ95_OBSBIASCOVARIANCE_H_
