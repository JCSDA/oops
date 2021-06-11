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

#include <memory>
#include <ostream>
#include <string>

#include <boost/noncopyable.hpp>

#include "lorenz95/ObservationTLAD.h"
#include "lorenz95/ObsTable.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace lorenz95 {
  class GomL95;
  class LocsL95;
  class ObsBias;
  class ObsDiags1D;
  class ObsVec1D;

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

  ObservationL95(const ObsTable &, const eckit::Configuration &);
  ~ObservationL95();

// Obs Operators
  void simulateObs(const GomL95 &, ObsVec1D &, const ObsBias &, ObsDiags1D &) const;

// Other
  const oops::Variables & requiredVars() const {return inputs_;}
  std::unique_ptr<LocsL95> locations() const;

  const ObsTable & table() const {return obsdb_;}

 private:
  void print(std::ostream &) const;
  const ObsTable & obsdb_;
  const oops::Variables inputs_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_OBSERVATIONL95_H_
