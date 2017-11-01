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

#include <ostream>
#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "model/ObsSpaceQG.h"
#include "model/ObservationsQG.h"
#include "model/ObsStreamTLAD.h"
#include "util/ObjectCounter.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
}

namespace qg {
  class GomQG;
  class ObsBias;
  class ObsBiasIncrement;
  class ObsVecQG;

// -----------------------------------------------------------------------------
/// Streamfunction observation for QG model.
/*!
 *  ObsStreamQG for QG model inherits from ObsEquivalent.
 */

class ObsStreamQG : public ObservationsQG,
                    private util::ObjectCounter<ObsStreamQG> {
 public:
  static const std::string classname() {return "qg::ObsStreamQG";}

  ObsStreamQG(ObsSpaceQG &, const eckit::Configuration &);
  virtual ~ObsStreamQG();

// Obs Operator
  void obsEquiv(const GomQG &, ObsVecQG &, const ObsBias &) const;

// Is there a way to put this in the TLAD class?
  LinearObsOp * getTLAD() const {return new ObsStreamTLAD(obsdb_, keyOperStrm_);}

// Other
  void generateObsError(const eckit::Configuration &);
  boost::shared_ptr<const VariablesQG> variables() const {return varin_;}

  int & toFortran() {return keyOperStrm_;}
  const int & toFortran() const {return keyOperStrm_;}

 private:
  void print(std::ostream &) const;
  ObsSpaceQG & obsdb_;
  const std::string obsname_;
  F90hop keyOperStrm_;
  boost::shared_ptr<const VariablesQG> varin_;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_OBSSTREAMQG_H_
