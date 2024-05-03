/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_OBSOPERATORTLAD_H_
#define QG_MODEL_OBSOPERATORTLAD_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "oops/qg/ObsDataQG.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace qg {
  class GomQG;
  class ObsBias;
  class ObsBiasIncrement;
  class ObsOpBaseTLAD;
  class ObsSpaceQG;
  class ObsVecQG;

// -----------------------------------------------------------------------------

class ObsOperatorTLAD : public util::Printable,
                        private boost::noncopyable {
 public:
  typedef ObsDataQG<int> QCFlags_;


  ObsOperatorTLAD(const ObsSpaceQG &, const eckit::Configuration &);

  ~ObsOperatorTLAD();

/// Obs Operator
  void setTrajectory(const GomQG &, const ObsBias &,
                     const QCFlags_ &);
  void simulateObsTL(const GomQG &, ObsVecQG &, const ObsBiasIncrement &,
                     const QCFlags_ &) const;
  void simulateObsAD(GomQG &, const ObsVecQG &, ObsBiasIncrement &,
                     const QCFlags_ &) const;

/// Other
  const oops::Variables & requiredVars() const;  // Required inputs requiredVars from Model

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ObsOpBaseTLAD> oper_;
};

// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_OBSOPERATORTLAD_H_
