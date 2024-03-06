/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_OBSWINDQG_H_
#define QG_MODEL_OBSWINDQG_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "oops/qg/ObsOpBaseQG.h"
#include "oops/qg/ObsSpaceQG.h"
#include "oops/qg/QgTraits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace qg {
  class GomQG;
  class LocationsQG;
  class ObsBias;
  class ObsVecQG;

// -----------------------------------------------------------------------------
/// Wind observation for QG model.

class ObsWindQG : public ObsOpBaseQG,
                  private util::ObjectCounter<ObsWindQG> {
 public:
  typedef ObsDataQG<int> QCFlags_;

  static const std::string classname() {return "qg::ObsWindQG";}

  ObsWindQG(const ObsSpaceQG &, const eckit::Configuration &);

// Obs Operators
  void simulateObs(const GomQG &, ObsVecQG &, const ObsBias &, const QCFlags_ &) const override;

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
#endif  // QG_MODEL_OBSWINDQG_H_
