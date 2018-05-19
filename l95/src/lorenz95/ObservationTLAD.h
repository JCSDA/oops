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

#include "lorenz95/L95Traits.h"
#include "oops/base/Variables.h"
#include "oops/interface/LinearObsOperBase.h"
#include "oops/util/ObjectCounter.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace lorenz95 {
  class GomL95;
  class ObsBias;
  class ObsBiasCorrection;
  class ObsTable;
  class ObsVec1D;

/// Observation for Lorenz 95 model.
/*!
 *  ObservationTLAD defines the linearized ObsOperator for Lorenz 95 model.
 */

// -----------------------------------------------------------------------------

class ObservationTLAD : public oops::LinearObsOperBase<L95Traits>,
                        private util::ObjectCounter<ObservationTLAD> {
 public:
  static const std::string classname() {return "lorenz95::ObservationTLAD";}

  ObservationTLAD(const ObsTable &, const eckit::Configuration &);
  ~ObservationTLAD();

// Obs Operators
  void setTrajectory(const GomL95 &, const ObsBias &);
  void obsEquivTL(const GomL95 &, ObsVec1D &, const ObsBiasCorrection &) const;
  void obsEquivAD(GomL95 &, const ObsVec1D &, ObsBiasCorrection &) const;

// Other
  const oops::Variables & variables() const {return inputs_;}

 private:
  void print(std::ostream &) const;
  const oops::Variables inputs_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_OBSERVATIONTLAD_H_
