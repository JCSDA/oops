/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_LINEAROBSOP_H_
#define QG_MODEL_LINEAROBSOP_H_

#include <map>
#include <string>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/interface/ObsOperatorBase.h"

namespace qg {
  class GomQG;
  class ObsVecQG;
  class ObsBias;
  class ObsBiasIncrement;
  class VariablesQG;

// -----------------------------------------------------------------------------

class LinearObsOp : private boost::noncopyable {
 public:
  LinearObsOp() {}
  virtual ~LinearObsOp() {}

// Obs Operators
  virtual void setTrajectory(const GomQG &, const ObsBias &) =0;
  virtual void obsEquivTL(const GomQG &, ObsVecQG &, const ObsBiasIncrement &) const =0;
  virtual void obsEquivAD(GomQG &, const ObsVecQG &, ObsBiasIncrement &) const =0;

// Other
  virtual boost::shared_ptr<const VariablesQG> variables() const =0;

  virtual int & toFortran() =0;
  virtual const int & toFortran() const =0;
};

// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_LINEAROBSOP_H_
