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

#include <boost/shared_ptr.hpp>

#include "oops/interface/ObsOperatorBase.h"
#include "util/ObjectCounter.h"
#include "lorenz95/ObservationTLAD.h"
#include "lorenz95/L95Traits.h"

// Forward declarations
namespace eckit {
  class Configuration;
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

class ObservationL95 : public oops::ObsOperatorBase<L95Traits>,
                       private util::ObjectCounter<ObservationL95> {
 public:
  static const std::string classname() {return "lorenz95::ObservationL95";}

  ObservationL95(ObsTable &, const eckit::Configuration &);
  ~ObservationL95();

// Obs Operators
  void obsEquiv(const GomL95 &, ObsVec1D &, const ObsBias &) const;

// Other
  void generateObsError(const eckit::Configuration &);
  boost::shared_ptr<const NoVariables> variables() const {return inputs_;}
  ObservationTLAD * newTLAD() const {return new ObservationTLAD(obsdb_);}

  const ObsTable & table() const {return obsdb_;}

 private:
  void print(std::ostream &) const;
  ObsTable & obsdb_;
  boost::shared_ptr<const NoVariables> inputs_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_OBSERVATIONL95_H_
