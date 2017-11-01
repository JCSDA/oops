/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_OBSERVATIONL95_H_
#define LORENZ95_OBSERVATIONL95_H_

#include <ostream>
#include <string>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "util/ObjectCounter.h"
#include "util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
}

namespace lorenz95 {
  class GomL95;
  class ObsBias;
  class ObsTable;
  class ObsVec1D;
  class NoVariables;

/// Observation for Lorenz 95 model.
/*!
 *  ObservationL95 defines ObsOperator for Lorenz 95 model.
 */

// -----------------------------------------------------------------------------

class ObservationL95 : public util::Printable,
                       private boost::noncopyable,
                       private util::ObjectCounter<ObservationL95> {
 public:
  static const std::string classname() {return "lorenz95::ObservationL95";}

  static ObservationL95 * create(ObsTable & ot, const eckit::Configuration & conf)
    {return new ObservationL95(ot, conf);}

  ~ObservationL95();

// Obs Operators
  void obsEquiv(const GomL95 &, ObsVec1D &, const ObsBias &) const;

// Other
  void generateObsError(const eckit::Configuration &);
  boost::shared_ptr<const NoVariables> variables() const {return inputs_;}

  const ObsTable & table() const {return obsdb_;}

 private:
  void print(std::ostream &) const;
  ObservationL95(ObsTable &, const eckit::Configuration &);
  ObsTable & obsdb_;
  boost::shared_ptr<const NoVariables> inputs_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_OBSERVATIONL95_H_
