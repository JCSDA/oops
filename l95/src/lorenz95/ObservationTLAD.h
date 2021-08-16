/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_OBSERVATIONTLAD_H_
#define LORENZ95_OBSERVATIONTLAD_H_

#include <ostream>
#include <string>

#include <boost/noncopyable.hpp>

#include "lorenz95/ObsTable.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace lorenz95 {
  class GomL95;
  class ObsBias;
  class ObsBiasCorrection;
  class ObsVec1D;

// -----------------------------------------------------------------------------

/// (Empty) parameters controlling the observation operator for the Lorenz 95 model.
class ObservationL95Parameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObservationL95Parameters, Parameters);
  // No parameters needed.
};

// -----------------------------------------------------------------------------

/// Observation for Lorenz 95 model.
/*!
 *  ObservationTLAD defines the linearized ObsOperator for Lorenz 95 model.
 */

class ObservationTLAD : public util::Printable,
                        private boost::noncopyable,
                        private util::ObjectCounter<ObservationTLAD> {
 public:
  typedef ObservationL95Parameters Parameters_;

  static const std::string classname() {return "lorenz95::ObservationTLAD";}

  ObservationTLAD(const ObsTable &, const Parameters_ &);

// Obs Operators
  void setTrajectory(const GomL95 &, const ObsBias &);
  void simulateObsTL(const GomL95 &, ObsVec1D &, const ObsBiasCorrection &) const;
  void simulateObsAD(GomL95 &, const ObsVec1D &, ObsBiasCorrection &) const;

// Other
  const oops::Variables & requiredVars() const {return inputs_;}

 private:
  void print(std::ostream &) const;
  const oops::Variables inputs_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_OBSERVATIONTLAD_H_
