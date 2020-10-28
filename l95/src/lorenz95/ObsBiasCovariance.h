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

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace lorenz95 {
  class ObsBias;
  class ObsBiasCorrection;
  class ObsTableView;

// -----------------------------------------------------------------------------

class ObsBiasCovariance : public util::Printable,
                          private boost::noncopyable,
                          private util::ObjectCounter<ObsBiasCovariance> {
 public:
  static const std::string classname() {return "lorenz95::ObsBiasCovariance";}

/// Constructor, destructor
  ObsBiasCovariance(const ObsTableView &, const eckit::Configuration &);
  ~ObsBiasCovariance() {}

/// Linear algebra operators
  void linearize(const ObsBias &, const eckit::Configuration &) {}
  void multiply(const ObsBiasCorrection &, ObsBiasCorrection &) const;
  void inverseMultiply(const ObsBiasCorrection &, ObsBiasCorrection &) const;
  void randomize(ObsBiasCorrection &) const;

  const eckit::Configuration & config() const {return conf_;}
  bool active() const {return active_;}

 private:
  void print(std::ostream &) const;
  const eckit::LocalConfiguration conf_;
  double variance_;
  bool active_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_OBSBIASCOVARIANCE_H_
