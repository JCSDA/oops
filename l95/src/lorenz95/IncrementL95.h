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

#ifndef LORENZ95_INCREMENTL95_H_
#define LORENZ95_INCREMENTL95_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "lorenz95/FieldL95.h"
#include "lorenz95/Resolution.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class LocalIncrement;
  class Variables;
}

namespace lorenz95 {
  class Iterator;
  class StateL95;

// -----------------------------------------------------------------------------

/// Increment Class: Difference between two states
/*!
 *  Some fields that are present in a State may not be present in
 *  an Increment. The Increment contains everything that is needed by
 *  the tangent-linear and adjoint models.
 */

// -----------------------------------------------------------------------------
class IncrementL95 : public util::Printable,
                     public util::Serializable,
                     private util::ObjectCounter<IncrementL95> {
 public:
  static const std::string classname() {return "lorenz95::IncrementL95";}

/// Constructor, destructor
  IncrementL95(const Resolution &, const oops::Variables &, const util::DateTime &);
  IncrementL95(const Resolution &, const IncrementL95 &);
  IncrementL95(const IncrementL95 &, const bool);
  virtual ~IncrementL95();

/// Basic operators
  void diff(const StateL95 &, const StateL95 &);
  void zero();
  void zero(const util::DateTime &);
  void ones();
  void dirac(const eckit::Configuration &);
  IncrementL95 & operator =(const IncrementL95 &);
  IncrementL95 & operator+=(const IncrementL95 &);
  IncrementL95 & operator-=(const IncrementL95 &);
  IncrementL95 & operator*=(const double &);
  void axpy(const double &, const IncrementL95 &, const bool check = true);
  double dot_product_with(const IncrementL95 &) const;
  void schur_product_with(const IncrementL95 &);
  void random(const size_t & seed = 1);

/// ATLAS
  void toFieldSet(atlas::FieldSet &) const {
    throw eckit::NotImplemented("IncrementL95::toFieldSet not implemented", Here());
  }
  void fromFieldSet(const atlas::FieldSet &) {
    throw eckit::NotImplemented("IncrementL95::fromFieldSet not implemented", Here());
  }

// Utilities
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm () const {return fld_.rms();}
  const util::DateTime & validTime() const {return time_;}
  util::DateTime & validTime() {return time_;}
  void updateTime(const util::Duration & dt) {time_ += dt;}

  oops::LocalIncrement getLocal(const Iterator &) const;
  void setLocal(const oops::LocalIncrement &, const Iterator &);

/// Access to data
  const FieldL95 & getField() const {return fld_;}
  FieldL95 & getField() {return fld_;}
  std::shared_ptr<const Resolution> geometry() const {
    std::shared_ptr<const Resolution> geom(new Resolution(fld_.resol()));
    return geom;
  }
  std::vector<double> & asVector() {return fld_.asVector();}
  const std::vector<double> & asVector() const {return fld_.asVector();}
  const oops::Variables & variables() const {return vars_;}

  void accumul(const double &, const StateL95 &);

/// Serialize and deserialize
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

 private:
  void print(std::ostream &) const override;
  FieldL95 fld_;
  util::DateTime time_;
  oops::Variables vars_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_INCREMENTL95_H_
