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

#include "lorenz95/L95Traits.h"
#include "lorenz95/ObservationTLAD.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsOperatorBase.h"
#include "util/ObjectCounter.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace lorenz95 {
  class GomL95;
  class ObsBias;
  class ObsTable;
  class ObsVec1D;

/// Observation for Lorenz 95 model.
/*!
 *  ObservationL95 defines ObsOperator for Lorenz 95 model.
 */

// -----------------------------------------------------------------------------

class ObservationL95 : public oops::ObsOperatorBase<L95Traits>,
                       private util::ObjectCounter<ObservationL95> {
 public:
  static const std::string classname() {return "lorenz95::ObservationL95";}

  ObservationL95(const ObsTable &, const eckit::Configuration &);
  ~ObservationL95();

// Obs Operators
  void obsEquiv(const GomL95 &, ObsVec1D &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const {return inputs_;}

  const ObsTable & table() const {return obsdb_;}

 private:
  void print(std::ostream &) const;
  const ObsTable & obsdb_;
  const oops::Variables inputs_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_OBSERVATIONL95_H_
