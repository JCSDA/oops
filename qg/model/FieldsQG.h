/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_FIELDSQG_H_
#define QG_MODEL_FIELDSQG_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/LocalIncrement.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

#include "oops/qg/GeometryQG.h"
#include "oops/qg/GeometryQGIterator.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class LocalIncrement;
}

namespace qg {
  class LocationsQG;
  class GomQG;

// -----------------------------------------------------------------------------
/// Class to represent a Fields for the QG model
class FieldsQG : public util::Printable,
                 public util::Serializable,
                 private util::ObjectCounter<FieldsQG> {
 public:
  static const std::string classname() {return "qg::FieldsQG";}

// Constructors and basic operators
  FieldsQG(const GeometryQG &, const oops::Variables &, const bool &, const util::DateTime &);
  FieldsQG(const FieldsQG &, const GeometryQG &, const bool ad = false);
  FieldsQG(const FieldsQG &, const oops::Variables &);
  FieldsQG(const FieldsQG &, const bool);
  FieldsQG(const FieldsQG &);
  ~FieldsQG();

  void zero();
  void zero(const util::DateTime &);
  void ones();
  FieldsQG & operator=(const FieldsQG &);
  FieldsQG & operator+=(const FieldsQG &);
  FieldsQG & operator-=(const FieldsQG &);
  FieldsQG & operator*=(const double &);
  void axpy(const double &, const FieldsQG &);
  double dot_product_with(const FieldsQG &) const;
  void schur_product_with(const FieldsQG &);
  void dirac(const eckit::Configuration &);
  void random();

// Interpolate full fields
  void changeResolution(const FieldsQG &);
  void add(const FieldsQG &);
  void diff(const FieldsQG &, const FieldsQG &);

// ATLAS FieldSet
  void toFieldSet(atlas::FieldSet &) const;
  void fromFieldSet(const atlas::FieldSet &);

// Utilities
  void read(const eckit::Configuration &);
  void analytic_init(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;
  std::shared_ptr<const GeometryQG> geometry() const {return geom_;}
  const oops::Variables & variables() const {return vars_;}
  oops::Variables & variables() {return vars_;}

  const util::DateTime & time() const {return time_;}
  util::DateTime & time() {return time_;}
  void updateTime(const util::Duration & dt) {time_ += dt;}

  const int & toFortran() const {return keyFlds_;}

  bool isForModel(const bool &) const;

  oops::LocalIncrement getLocal(const GeometryQGIterator &) const;
  void setLocal(const oops::LocalIncrement &, const GeometryQGIterator &);

/// Serialization
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

 private:
  void print(std::ostream &) const override;
  F90flds keyFlds_;
  std::shared_ptr<const GeometryQG> geom_;
  oops::Variables vars_;
  const bool lbc_;
  util::DateTime time_;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_FIELDSQG_H_
