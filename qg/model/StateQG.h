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

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "model/FieldsQG.h"
#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace qg {
  class GomQG;
  class LocationsQG;
  class GeometryQG;
  class IncrementQG;
  class VariablesQG;

/// QG model state
/*!
 * A State contains everything that is needed to propagate the state
 * forward in time.
 */

// -----------------------------------------------------------------------------
class StateQG : public util::Printable,
                private util::ObjectCounter<StateQG> {
 public:
  static const std::string classname() {return "qg::StateQG";}

/// Constructor, destructor
  StateQG(const GeometryQG &, const VariablesQG &, const util::DateTime &);  // Is it used?
  StateQG(const GeometryQG &, const eckit::Configuration &);
  StateQG(const GeometryQG &, const StateQG &);
  StateQG(const StateQG &);
  virtual ~StateQG();
  StateQG & operator=(const StateQG &);

/// Interpolate to observation location
  void interpolate(const LocationsQG &, const VariablesQG &, GomQG &) const;

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
  void activateModel();
  void deactivateModel();

  void zero();
  void accumul(const double &, const StateQG &);

 private:
  void print(std::ostream &) const;
  boost::scoped_ptr<FieldsQG> fields_;
  boost::scoped_ptr<FieldsQG> stash_;
};
// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_STATEQG_H_
