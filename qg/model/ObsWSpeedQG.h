/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_OBSWSPEEDQG_H_
#define QG_MODEL_OBSWSPEEDQG_H_

#include <ostream>
#include <string>

#include <boost/shared_ptr.hpp>

#include "model/ObsSpaceQG.h"
#include "model/QgTraits.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsOperatorBase.h"
#include "oops/util/ObjectCounter.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace qg {
  class GomQG;
  class LocationsQG;
  class ObsBias;
  class ObsBiasIncrement;
  class ObsVecQG;

// -----------------------------------------------------------------------------
/// Wind speed observation for QG model.

class ObsWSpeedQG : public oops::ObsOperatorBase<QgTraits>,
                    private util::ObjectCounter<ObsWSpeedQG> {
 public:
  static const std::string classname() {return "qg::ObsWSpeedQG";}

  ObsWSpeedQG(const ObsSpaceQG &, const eckit::Configuration &);
  virtual ~ObsWSpeedQG();

// Obs Operator
  void obsEquiv(const GomQG &, ObsVecQG &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const override {return varin_;}

  int & toFortran() {return keyOperWspeed_;}
  const int & toFortran() const {return keyOperWspeed_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperWspeed_;
  const oops::Variables varin_;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_OBSWSPEEDQG_H_
