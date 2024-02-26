/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_OBSBIASINCREMENT_H_
#define QG_MODEL_OBSBIASINCREMENT_H_

#include <iostream>
#include <vector>

#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

namespace eckit {
  class Configuration;
}

namespace qg {
  class ObsBias;
  class ObsSpaceQG;

// -----------------------------------------------------------------------------

class ObsBiasIncrement : public util::Printable,
                         public util::Serializable {
 public:
/// Constructor, destructor
  ObsBiasIncrement();
  ObsBiasIncrement(const ObsSpaceQG &, const eckit::Configuration &);
  ObsBiasIncrement(const ObsBiasIncrement &, const bool copy = true);

/// Linear algebra operators
  void diff(const ObsBias &, const ObsBias &);
  void zero();
  ObsBiasIncrement & operator=(const ObsBiasIncrement &);
  ObsBiasIncrement & operator+=(const ObsBiasIncrement &);
  ObsBiasIncrement & operator-=(const ObsBiasIncrement &);
  ObsBiasIncrement & operator*=(const double);
  void axpy(const double, const ObsBiasIncrement &);
  double dot_product_with(const ObsBiasIncrement &) const;

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const;

  double & operator[](const unsigned int ii) {return bias_[ii];}
  const double & operator[](const unsigned int ii) const {return bias_[ii];}

  double & stream() {return bias_[0];}
  double & wind() {return bias_[1];}
  double & wspd() {return bias_[3];}
  const double & stream() const {return bias_[0];}
  const double & wind() const {return bias_[1];}
  const double & wspd() const {return bias_[3];}

/// Serialization
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

 private:
  void print(std::ostream &) const override;
  void makePassive();

  std::vector<double> bias_;
  std::vector<bool> active_;
};

// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_OBSBIASINCREMENT_H_
