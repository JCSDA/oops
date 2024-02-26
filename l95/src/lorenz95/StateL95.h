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

#ifndef LORENZ95_STATEL95_H_
#define LORENZ95_STATEL95_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "lorenz95/FieldL95.h"
#include "lorenz95/Resolution.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace lorenz95 {
  class IncrementL95;

/// L95 model state
/*!
 * A State contains everything that is needed to propagate the state
 * forward in time.
 */

// -----------------------------------------------------------------------------
class StateL95 : public util::Printable,
                 private util::ObjectCounter<StateL95> {
 public:
  static const std::string classname() {return "lorenz95::StateL95";}

/// Constructor, destructor
  StateL95(const Resolution &, const oops::Variables &, const util::DateTime &);
  StateL95(const Resolution &, const eckit::Configuration &);
  StateL95(const Resolution &, const StateL95 &);
  StateL95(const oops::Variables &, const StateL95 &);
  StateL95(const StateL95 &);
  virtual ~StateL95();
  StateL95 & operator=(const StateL95 &);

/// Interactions with increments
  StateL95 & operator+=(const IncrementL95 &);

// Utilities
  const FieldL95 & getField() const {return fld_;}
  FieldL95 & getField() {return fld_;}
  std::shared_ptr<const Resolution> geometry() const {
    std::shared_ptr<const Resolution> geom(new Resolution(fld_.resol()));
    return geom;
  }

/// ATLAS
  void toFieldSet(atlas::FieldSet &) const {
    throw eckit::NotImplemented("StateL95::toFieldSet not implemented", Here());
  }
  void fromFieldSet(const atlas::FieldSet &) {
    throw eckit::NotImplemented("StateL95::fromFieldSet not implemented", Here());
  }

  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm () const {return fld_.rms();}
  const util::DateTime & validTime() const {return time_;}
  void updateTime(const util::Duration & dt) {time_ += dt;}
  util::DateTime & validTime() {return time_;}
  const oops::Variables & variables() const {return vars_;}

// For accumulator
  void zero();
  void accumul(const double &, const StateL95 &);

/// Serialize and deserialize
  size_t serialSize() const;
  void serialize(std::vector<double> &) const;
  void deserialize(const std::vector<double> &, size_t &);

 private:
  void print(std::ostream &) const;

  FieldL95 fld_;
  util::DateTime time_;
  oops::Variables vars_;
};
// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_STATEL95_H_
