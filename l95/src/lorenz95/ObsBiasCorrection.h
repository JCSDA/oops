/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_OBSBIASCORRECTION_H_
#define LORENZ95_OBSBIASCORRECTION_H_

#include <cmath>
#include <iostream>

#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace lorenz95 {
  class ObsBias;

// -----------------------------------------------------------------------------

class ObsBiasCorrection : public util::Printable {
 public:
/// Constructor, destructor
  explicit ObsBiasCorrection();
  explicit ObsBiasCorrection(const eckit::Configuration &);
  ObsBiasCorrection(const ObsBiasCorrection &, const bool copy = true);
  ObsBiasCorrection(const ObsBiasCorrection &, const eckit::Configuration &);
  ~ObsBiasCorrection() {}

/// Linear algebra operators
  void diff(const ObsBias &, const ObsBias &);
  void zero();
  ObsBiasCorrection & operator=(const ObsBiasCorrection &);
  ObsBiasCorrection & operator+=(const ObsBiasCorrection &);
  ObsBiasCorrection & operator-=(const ObsBiasCorrection &);
  ObsBiasCorrection & operator*=(const double);
  void axpy(const double, const ObsBiasCorrection &);
  double dot_product_with(const ObsBiasCorrection &) const;

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const {return std::abs(bias_);}

  double & value() {return bias_;}
  const double & value() const {return bias_;}

 private:
  void print(std::ostream &) const;
  double bias_;
  bool active_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_OBSBIASCORRECTION_H_
