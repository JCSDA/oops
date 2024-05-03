/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_OBSSTREAMTLAD_H_
#define QG_MODEL_OBSSTREAMTLAD_H_

#include <string>

#include "oops/base/Variables.h"
#include "oops/qg/ObsOpBaseTLAD.h"
#include "oops/util/ObjectCounter.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace qg {
  class GomQG;
  class ObsBias;
  class ObsBiasIncrement;
  class ObsSpaceQG;
  class ObsVecQG;

// -----------------------------------------------------------------------------
/// Streamfunction TL/AD observation operator for QG model.

class ObsStreamTLAD : public ObsOpBaseTLAD,
                      private util::ObjectCounter<ObsStreamTLAD> {
 public:
  typedef ObsDataQG<int> QCFlags_;
  static const std::string classname() {return "qg::ObsStreamTLAD";}

  ObsStreamTLAD(const ObsSpaceQG &, const eckit::Configuration &);

// Obs Operators
  void setTrajectory(const GomQG &, const ObsBias &,
                     const QCFlags_ &) override;
  void simulateObsTL(const GomQG &, ObsVecQG &, const ObsBiasIncrement &,
                     const QCFlags_ &) const override;
  void simulateObsAD(GomQG &, const ObsVecQG &, ObsBiasIncrement &,
                     const QCFlags_ &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

 private:
  void print(std::ostream &) const override;
  const ObsSpaceQG & obsdb_;
  const oops::Variables varin_;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_OBSSTREAMTLAD_H_
