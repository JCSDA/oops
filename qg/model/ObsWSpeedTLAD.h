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
#include <vector>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "model/GomQG.h"
#include "model/LinearObsOp.h"
#include "model/ObsSpaceQG.h"
#include "util/ObjectCounter.h"

// Forward declarations
namespace util {
  class DateTime;
}

namespace qg {
  class ObsBias;
  class ObsBiasIncrement;
  class ObsVecQG;

// -----------------------------------------------------------------------------

class ObsWSpeedTLAD : public LinearObsOp, private util::ObjectCounter<ObsWSpeedTLAD> {
 public:
  static const std::string classname() {return "qg::ObsWSpeedTLAD";}

  ObsWSpeedTLAD(const ObsSpaceQG &, const int &);
  virtual ~ObsWSpeedTLAD();

// Obs Operators
  void setTrajectory(const GomQG &, const ObsBias &);
  void obsEquivTL(const GomQG &, ObsVecQG &, const ObsBiasIncrement &) const;
  void obsEquivAD(GomQG &, const ObsVecQG &, ObsBiasIncrement &) const;

// Other
  boost::shared_ptr<const VariablesQG> variables() const {return varin_;}

  int & toFortran() {return keyOperWspeed_;}
  const int & toFortran() const {return keyOperWspeed_;}

 private:
  F90hop keyOperWspeed_;
  GomQG traj_;
  boost::shared_ptr<const VariablesQG> varin_;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_OBSWSPEEDTLAD_H_
