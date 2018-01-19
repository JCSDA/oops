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

#include <ostream>
#include <string>

#include <boost/shared_ptr.hpp>

#include "model/GeometryQG.h"
#include "model/VariablesQG.h"
#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class UnstructuredGrid;
  class Variables;
}

namespace qg {
  class LocationsQG;
  class GomQG;

// -----------------------------------------------------------------------------
/// Class to represent a FieldSet for the QG model
class FieldsQG : public util::Printable,
                 private util::ObjectCounter<FieldsQG> {
 public:
  static const std::string classname() {return "qg::FieldsQG";}

// Constructors and basic operators
  FieldsQG(const GeometryQG &, const oops::Variables &, const util::DateTime &);
  FieldsQG(const FieldsQG &, const GeometryQG &);
  FieldsQG(const FieldsQG &, const oops::Variables &);
  FieldsQG(const FieldsQG &, const bool);
  FieldsQG(const FieldsQG &);
  ~FieldsQG();

  void zero();
  void zero(const util::DateTime &);
  FieldsQG & operator=(const FieldsQG &);
  FieldsQG & operator+=(const FieldsQG &);
  FieldsQG & operator-=(const FieldsQG &);
  FieldsQG & operator*=(const double &);
  void axpy(const double &, const FieldsQG &);
  double dot_product_with(const FieldsQG &) const;
  void schur_product_with(const FieldsQG &);
  void dirac(const eckit::Configuration &);
  void random();

// Interpolate to given location
  void interpolate(const LocationsQG &, const oops::Variables &, GomQG &) const;
  void interpolateTL(const LocationsQG &, const oops::Variables &, GomQG &) const;
  void interpolateAD(const LocationsQG &, const oops::Variables &, const GomQG &);

// Interpolate full fields
  void changeResolution(const FieldsQG &);
  void add(const FieldsQG &);
  void diff(const FieldsQG &, const FieldsQG &);

// Define and convert to/from unstructured grid
  void define(oops::UnstructuredGrid &) const;
  void convert_to(oops::UnstructuredGrid &) const;
  void convert_from(const oops::UnstructuredGrid &);

// Utilities
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;
  boost::shared_ptr<const GeometryQG> geometry() const {return geom_;}

  const util::DateTime & time() const {return time_;}
  util::DateTime & time() {return time_;}

  int & toFortran() {return keyFlds_;}
  const int & toFortran() const {return keyFlds_;}

  bool isForModel(const bool) const;

 private:
  void print(std::ostream &) const;
  F90flds keyFlds_;
  boost::shared_ptr<const GeometryQG> geom_;
  const VariablesQG vars_;
  util::DateTime time_;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_FIELDSQG_H_
