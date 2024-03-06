/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_OBSSTREAMQG_H_
#define QG_MODEL_OBSSTREAMQG_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "oops/qg/ObsDataQG.h"
#include "oops/qg/ObsOpBaseQG.h"
#include "oops/qg/ObsSpaceQG.h"
#include "oops/qg/QgTraits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace qg {
  class GomQG;
  class ObsBias;
  class ObsVecQG;

// -----------------------------------------------------------------------------
/// Streamfunction observation for QG model.

class ObsStreamQG : public ObsOpBaseQG,
                    private util::ObjectCounter<ObsStreamQG> {
 public:
  static const std::string classname() {return "qg::ObsStreamQG";}
  typedef ObsDataQG<int> QCFlags_;
  ObsStreamQG(const ObsSpaceQG &, const eckit::Configuration &);

// Obs Operator
  void simulateObs(const GomQG &, ObsVecQG &, const ObsBias &,
                   const QCFlags_ &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}
  Locations_ locations() const override;

 private:
  void print(std::ostream &) const override;
  const ObsSpaceQG & obsdb_;
  const oops::Variables varin_;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_OBSSTREAMQG_H_
