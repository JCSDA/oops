/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_OBSOPERATORQG_H_
#define QG_MODEL_OBSOPERATORQG_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "model/ObsOperatorParameters.h"

#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace oops {
  template <typename OBS> class Locations;
}

namespace qg {
  class GomQG;
  class ObsBias;
  class ObsDiagsQG;
  class ObsOpBaseQG;
  class ObsSpaceQG;
  class ObsVecQG;
  struct QgObsTraits;

// -----------------------------------------------------------------------------

class ObsOperatorQG : public util::Printable,
                      private boost::noncopyable {
 public:
  typedef ObservationParameters Parameters_;
  typedef oops::Locations<QgObsTraits> Locations_;

  ObsOperatorQG(const ObsSpaceQG &, const Parameters_ &);
  ~ObsOperatorQG();

/// Obs Operator
  void simulateObs(const GomQG &, ObsVecQG &, const ObsBias &, ObsVecQG &, ObsDiagsQG &) const;

/// Other
  const oops::Variables & requiredVars() const;  // Required input requiredVars from Model
/// Model variable interpolation paths
  Locations_ locations() const;

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ObsOpBaseQG> oper_;
};

// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_OBSOPERATORQG_H_
