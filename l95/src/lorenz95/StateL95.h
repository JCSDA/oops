/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_STATEL95_H_
#define LORENZ95_STATEL95_H_

#include <ostream>
#include <string>

#include "lorenz95/FieldL95.h"
#include "lorenz95/Resolution.h"

#include "oops/base/GridPoint.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class GridPoint;
  class UnstructuredGrid;
  class Variables;
}

namespace lorenz95 {
  class GomL95;
  class IncrementL95;
  class Iterator;
  class LocsL95;
  class ModelBias;
  class ModelL95;
  class ModelTrajectory;
  class Nothing;

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
  StateL95(const Resolution &, const oops::Variables &, const eckit::Configuration &);
  StateL95(const Resolution &, const StateL95 &);
  StateL95(const StateL95 &);
  virtual ~StateL95();
  StateL95 & operator=(const StateL95 &);

/// Get state values at obs locations
  void getValues(const LocsL95 &, const oops::Variables &, GomL95 &) const;
  void getValues(const LocsL95 &, const oops::Variables &, GomL95 &, Nothing &) const;

/// Interactions with increments
  StateL95 & operator+=(const IncrementL95 &);

// Utilities
  const FieldL95 & getField() const {return fld_;}
  FieldL95 & getField() {return fld_;}
  boost::shared_ptr<const Resolution> geometry() const {
    boost::shared_ptr<const Resolution> geom(new Resolution(fld_.resol()));
    return geom;
  }

  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm () const {return fld_.rms();}
  const util::DateTime & validTime() const {return time_;}
  util::DateTime & validTime() {return time_;}

  oops::GridPoint getPoint(const Iterator &) const;

// For accumulator
  void zero();
  void accumul(const double &, const StateL95 &);

 private:
  void print(std::ostream &) const;
  FieldL95 fld_;
  util::DateTime time_;
};
// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_STATEL95_H_
