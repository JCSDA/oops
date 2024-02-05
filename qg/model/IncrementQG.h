/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_INCREMENTQG_H_
#define QG_MODEL_INCREMENTQG_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/LocalIncrement.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

#include "oops/qg/FieldsQG.h"
#include "oops/qg/GeometryQG.h"
#include "oops/qg/GeometryQGIterator.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class LocalIncrement;
  class Variables;
}

namespace qg {
  class GomQG;
  class LocationsQG;
  class GeometryQG;
  class ModelBiasIncrement;
  class ErrorCovarianceQG;
  class StateQG;

/// Increment Class: Difference between two states
/*!
 *  Some fields that are present in a State may not be present in
 *  an Increment. The Increment contains everything that is needed by
 *  the tangent-linear and adjoint models.
 */

// -----------------------------------------------------------------------------

class IncrementQG : public util::Printable,
                    public util::Serializable,
                    private util::ObjectCounter<IncrementQG> {
 public:
  static const std::string classname() {return "qg::IncrementQG";}

/// Constructor, destructor
  IncrementQG(const GeometryQG &, const oops::Variables &, const util::DateTime &);
  IncrementQG(const GeometryQG &, const IncrementQG &, const bool ad);
  IncrementQG(const IncrementQG &, const bool);
  IncrementQG(const IncrementQG &);
  virtual ~IncrementQG();

/// Basic operators
  void diff(const StateQG &, const StateQG &);
  void zero();
  void zero(const util::DateTime &);
  void ones();
  IncrementQG & operator =(const IncrementQG &);
  IncrementQG & operator+=(const IncrementQG &);
  IncrementQG & operator-=(const IncrementQG &);
  IncrementQG & operator*=(const double &);
  void axpy(const double &, const IncrementQG &, const bool check = true);
  double dot_product_with(const IncrementQG &) const;
  void schur_product_with(const IncrementQG &);
  void random();
  void dirac(const eckit::Configuration &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const {return fields_->norm();}
  const util::DateTime & validTime() const {return fields_->time();}
  util::DateTime & validTime() {return fields_->time();}
  void updateTime(const util::Duration & dt) {fields_->time() += dt;}

/// ATLAS FieldSet
  void toFieldSet(atlas::FieldSet &) const;
  void fromFieldSet(const atlas::FieldSet &);

/// Access to fields
  FieldsQG & fields() {return *fields_;}
  const FieldsQG & fields() const {return *fields_;}

  std::shared_ptr<const GeometryQG> geometry() const {
    return fields_->geometry();
  }

/// Other
  void accumul(const double &, const StateQG &);
  oops::LocalIncrement getLocal(const GeometryQGIterator &) const;
  void setLocal(const oops::LocalIncrement &, const GeometryQGIterator &);

/// Serialization
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

  const oops::Variables & variables() const {return fields_->variables();}

/// Data
 private:
  void print(std::ostream &) const override;
  const bool lbc_ = false;
  std::unique_ptr<FieldsQG> fields_;
};
// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_INCREMENTQG_H_
