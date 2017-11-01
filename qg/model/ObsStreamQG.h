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

#include <boost/shared_ptr.hpp>

#include "oops/interface/ObsOperatorBase.h"
#include "model/ObsSpaceQG.h"
#include "model/QgTraits.h"
#include "util/ObjectCounter.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace qg {
  class GomQG;
  class ObsBias;
  class ObsBiasIncrement;
  class ObsVecQG;

// -----------------------------------------------------------------------------
/// Streamfunction observation for QG model.

class ObsStreamQG : public oops::ObsOperatorBase<QgTraits>,
                    private util::ObjectCounter<ObsStreamQG> {
 public:
  static const std::string classname() {return "qg::ObsStreamQG";}

  ObsStreamQG(const ObsSpaceQG &, const eckit::Configuration &);
  virtual ~ObsStreamQG();

// Obs Operator
  void obsEquiv(const GomQG &, ObsVecQG &, const ObsBias &) const;

// Other
  boost::shared_ptr<const VariablesQG> variables() const {return varin_;}

  int & toFortran() {return keyOperStrm_;}
  const int & toFortran() const {return keyOperStrm_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperStrm_;
  boost::shared_ptr<const VariablesQG> varin_;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_OBSSTREAMQG_H_
