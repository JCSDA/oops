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
#include <vector>

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "oops/qg/FieldsQG.h"

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

/// QG model state
// -----------------------------------------------------------------------------
class StateQG : public util::Printable,
                private util::ObjectCounter<StateQG> {
 public:
  static const std::string classname() {return "qg::StateQG";}

/// Constructor, destructor
  StateQG(const GeometryQG &, const oops::Variables &, const util::DateTime &);  // Is it used?
  StateQG(const GeometryQG &, const eckit::Configuration &);
  StateQG(const GeometryQG &, const StateQG &);
  StateQG(const oops::Variables &, const StateQG &);
  StateQG(const StateQG &);
  virtual ~StateQG();
  StateQG & operator=(const StateQG &);

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
  void updateTime(const util::Duration & dt) {fields_->updateTime(dt);}

/// Access to fields
  FieldsQG & fields() {return *fields_;}
  const FieldsQG & fields() const {return *fields_;}
  std::shared_ptr<const GeometryQG> geometry() const {
    return fields_->geometry();
  }
  const oops::Variables & variables() const {return fields_->variables();}

/// Serialization
  size_t serialSize() const;
  void serialize(std::vector<double> &) const;
  void deserialize(const std::vector<double> &, size_t &);

/// ATLAS FieldSet
  void toFieldSet(atlas::FieldSet &) const;
  void fromFieldSet(const atlas::FieldSet &);

/// Other
  void zero();
  void accumul(const double &, const StateQG &);

 private:
  void print(std::ostream &) const;
  const bool lbc_ = true;
  std::unique_ptr<FieldsQG> fields_;
};
// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_STATEQG_H_
