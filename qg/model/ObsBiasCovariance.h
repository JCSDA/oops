/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_OBSBIASCOVARIANCE_H_
#define QG_MODEL_OBSBIASCOVARIANCE_H_

#include <array>
#include <memory>
#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "model/ObsBias.h"
#include "model/ObsBiasPreconditioner.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace qg {
  class ObsBias;
  class ObsBiasIncrement;
  class ObsSpaceQG;

// -----------------------------------------------------------------------------

class ObsBiasCovariance : public util::Printable,
                          private boost::noncopyable,
                          private util::ObjectCounter<ObsBiasCovariance> {
 public:
  static const std::string classname() {return "qg::ObsBiasCovariance";}

/// Constructor, destructor
  ObsBiasCovariance(const ObsSpaceQG &, const eckit::Configuration &);
  ~ObsBiasCovariance() {}

/// Linear algebra operators
  void linearize(const ObsBias &, const eckit::Configuration &) {}
  void multiply(const ObsBiasIncrement &, ObsBiasIncrement &) const;
  void inverseMultiply(const ObsBiasIncrement &, ObsBiasIncrement &) const;
  void randomize(ObsBiasIncrement &) const;

  std::unique_ptr<ObsBiasPreconditioner> preconditioner() const;

/// I/O and diagnostics
  void write(const eckit::Configuration &) const {}
  bool active(const unsigned int ii) const {return variance_[ii] > 0.0;}

 private:
  void print(std::ostream &) const;
  std::array<double, ObsBias::ntypes> variance_;
};

// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_OBSBIASCOVARIANCE_H_
