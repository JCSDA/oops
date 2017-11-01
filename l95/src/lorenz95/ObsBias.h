/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_OBSBIAS_H_
#define LORENZ95_OBSBIAS_H_

#include <cmath>
#include <iostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace lorenz95 {
  class ObsBiasCorrection;

/// Class to handle observation bias parameters.

// -----------------------------------------------------------------------------

class ObsBias : public util::Printable,
                private boost::noncopyable,
                private util::ObjectCounter<ObsBias> {
 public:
  static const std::string classname() {return "lorenz95::ObsBias";}

  explicit ObsBias(const eckit::Configuration &);
  ObsBias(const ObsBias &, const bool);
  ~ObsBias() {}

  ObsBias & operator+=(const ObsBiasCorrection &);

  const double & value() const {return bias_;}
  double & value() {return bias_;}

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const {return std::abs(bias_);}

 private:
  void print(std::ostream &) const;
  double bias_;
  bool active_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_OBSBIAS_H_
