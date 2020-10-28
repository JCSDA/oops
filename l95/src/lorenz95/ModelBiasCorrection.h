/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_MODELBIASCORRECTION_H_
#define LORENZ95_MODELBIASCORRECTION_H_

#include <cmath>
#include <iostream>
#include <vector>

#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

namespace eckit {
  class Configuration;
}

namespace lorenz95 {
  class ModelBias;
  class ModelBiasCovariance;
  class Resolution;

// -----------------------------------------------------------------------------

class ModelBiasCorrection : public util::Printable,
                            public util::Serializable {
 public:
/// Constructor, destructor
  ModelBiasCorrection(const Resolution &, const eckit::Configuration &);
  ModelBiasCorrection(const ModelBiasCorrection &, const bool);
  ModelBiasCorrection(const ModelBiasCorrection &, const eckit::Configuration &);
  ~ModelBiasCorrection() {}

/// Linear algebra operators
  void diff(const ModelBias &, const ModelBias &);
  void zero();
  ModelBiasCorrection & operator=(const ModelBiasCorrection &);
  ModelBiasCorrection & operator+=(const ModelBiasCorrection &);
  ModelBiasCorrection & operator-=(const ModelBiasCorrection &);
  ModelBiasCorrection & operator*=(const double);
  void axpy(const double, const ModelBiasCorrection &);
  double dot_product_with(const ModelBiasCorrection &) const;

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const {return std::abs(bias_);}

  double & bias() {return bias_;}
  const double & bias() const {return bias_;}

/// Serialize and deserialize
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

 private:
  ModelBiasCorrection(const ModelBiasCorrection &);
  void print(std::ostream &) const override;
  double bias_;
  bool active_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_MODELBIASCORRECTION_H_
