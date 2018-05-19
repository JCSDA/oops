/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_OBSWSPEEDTLAD_H_
#define QG_MODEL_OBSWSPEEDTLAD_H_

#include <string>

#include <boost/shared_ptr.hpp>

#include "model/QgTraits.h"
#include "oops/base/Variables.h"
#include "oops/interface/LinearObsOperBase.h"
#include "util/ObjectCounter.h"

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
/// Wind speed TL/AD observation operator for QG model.

class ObsWSpeedTLAD : public oops::LinearObsOperBase<QgTraits>,
                      private util::ObjectCounter<ObsWSpeedTLAD> {
 public:
  static const std::string classname() {return "qg::ObsWSpeedTLAD";}

  ObsWSpeedTLAD(const ObsSpaceQG &, const eckit::Configuration &);
  virtual ~ObsWSpeedTLAD();

// Obs Operators
  void setTrajectory(const GomQG &, const ObsBias &) override;
  void obsEquivTL(const GomQG &, ObsVecQG &, const ObsBiasIncrement &) const override;
  void obsEquivAD(GomQG &, const ObsVecQG &, ObsBiasIncrement &) const override;

// Other
  const oops::Variables & variables() const override {return varin_;}

  int & toFortran() {return keyOperWspeed_;}
  const int & toFortran() const {return keyOperWspeed_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperWspeed_;
  GomQG traj_;
  const oops::Variables varin_;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_OBSWSPEEDTLAD_H_
