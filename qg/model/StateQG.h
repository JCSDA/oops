/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_STATEQG_H_
#define QG_MODEL_STATEQG_H_

#include <memory>
#include <ostream>
#include <string>


#include "model/FieldsQG.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace qg {
  class GomQG;
  class LocationsQG;
  class GeometryQG;
  class IncrementQG;
  class Nothing;

/// QG model state
// -----------------------------------------------------------------------------
class StateQG : public util::Printable,
                private util::ObjectCounter<StateQG> {
 public:
  static const std::string classname() {return "qg::StateQG";}

/// Constructor, destructor
  StateQG(const GeometryQG &, const oops::Variables &, const util::DateTime &);  // Is it used?
  StateQG(const GeometryQG &, const oops::Variables &, const eckit::Configuration &);
  StateQG(const GeometryQG &, const StateQG &);
  StateQG(const StateQG &);
  virtual ~StateQG();
  StateQG & operator=(const StateQG &);

/// Get state values at observation locations
  void getValues(const LocationsQG &, const oops::Variables &, GomQG &) const;
  void getValues(const LocationsQG &, const oops::Variables &, GomQG &, Nothing &) const;

/// Interpolate full fields
  void changeResolution(const StateQG & xx);

/// Interactions with Increment
  StateQG & operator+=(const IncrementQG &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const {return fields_->norm();}
  const util::DateTime & validTime() const {return fields_->time();}
  util::DateTime & validTime() {return fields_->time();}

/// Access to fields
  FieldsQG & fields() {return *fields_;}
  const FieldsQG & fields() const {return *fields_;}
  boost::shared_ptr<const GeometryQG> geometry() const {
    return fields_->geometry();
  }

/// Other
  void zero();
  void accumul(const double &, const StateQG &);

 private:
  void print(std::ostream &) const;
  const bool lbc_ = true;
  std::unique_ptr<FieldsQG> fields_;
  std::unique_ptr<FieldsQG> stash_;
};
// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_STATEQG_H_
