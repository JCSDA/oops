/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_FIELDL95_H_
#define LORENZ95_FIELDL95_H_

#include <ostream>
#include <string>
#include <vector>

#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

namespace eckit {
  class Configuration;
}

namespace lorenz95 {
// Forward declarations
  class GomL95;
  class Resolution;

// -----------------------------------------------------------------------------
/// Class to represent fields for the L95 model
class FieldL95 : public util::Printable,
                 public util::Serializable {
 public:
  static const std::string classname() {return "lorenz95::FieldL95";}

/// Constructors and basic operators
  explicit FieldL95(const Resolution &);
  FieldL95(const FieldL95 &, const Resolution &);
  explicit FieldL95(const FieldL95 &, const bool copy = true);
  ~FieldL95() {}

/// Linear algebra
  void zero();
  void ones();
  void dirac(const eckit::Configuration &);
  FieldL95 & operator=(const FieldL95 &);
  FieldL95 & operator+=(const FieldL95 &);
  FieldL95 & operator-=(const FieldL95 &);
  FieldL95 & operator*=(const double &);
  void diff(const FieldL95 &, const FieldL95 &);
  void axpy(const double &, const FieldL95 &);
  double dot_product_with(const FieldL95 &) const;
  void schur(const FieldL95 &);
  void random(const size_t &);
  void generate(const eckit::Configuration &);

/// Utilities
  void read(std::ifstream &);
  void write(std::ofstream &) const;
  double rms() const;

/// Set and get
  const int & resol() const {return resol_;}
  double & operator[](const int ii) {return x_[ii];}
  const double & operator[](const int ii) const {return x_[ii];}
  std::vector<double> & asVector() {return x_;}
  const std::vector<double> & asVector() const {return x_;}

/// Serialize and deserialize
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

 private:
  void print(std::ostream &) const override;
  const int resol_;
  std::vector<double> x_;
};
// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_FIELDL95_H_
