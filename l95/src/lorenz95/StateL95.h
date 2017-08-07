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

#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class UnstructuredGrid;
}

namespace lorenz95 {
  class GomL95;
  class IncrementL95;
  class LocsL95;
  class ModelBias;
  class ModelL95;
  class ModelTrajectory;
  class NoVariables;

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
  StateL95(const Resolution &, const NoVariables &, const util::DateTime &);
  StateL95(const Resolution &, const eckit::Configuration &);
  StateL95(const Resolution &, const StateL95 &);
  StateL95(const StateL95 &);
  virtual ~StateL95();
  StateL95 & operator=(const StateL95 &);

/// Interpolate to observation location
  void interpolate(const LocsL95 &, GomL95 &) const;

/// Interactions with increments
  StateL95 & operator+=(const IncrementL95 &);

/// Convert to/from generic unstructured grid
  void convert_to(oops::UnstructuredGrid &) const;
  void convert_from(const oops::UnstructuredGrid &);

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

 private:
  void print(std::ostream &) const;
  FieldL95 fld_;
  util::DateTime time_;
};
// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_STATEL95_H_
