/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_OBSWINDTLAD_H_
#define QG_MODEL_OBSWINDTLAD_H_

#include <string>

#include <boost/shared_ptr.hpp>

#include "oops/interface/LinearObsOperBase.h"
#include "util/ObjectCounter.h"
#include "model/QgTraits.h"

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
/// Wind TL/AD observation operator for QG model.

class ObsWindTLAD : public oops::LinearObsOperBase<QgTraits>,
                    private util::ObjectCounter<ObsWindTLAD> {
 public:
  static const std::string classname() {return "qg::ObsWindTLAD";}

  ObsWindTLAD(const ObsSpaceQG &, const eckit::Configuration &);
  virtual ~ObsWindTLAD();

// Obs Operators
  void setTrajectory(const GomQG &, const ObsBias &);
  void obsEquivTL(const GomQG &, ObsVecQG &, const ObsBiasIncrement &) const;
  void obsEquivAD(GomQG &, const ObsVecQG &, ObsBiasIncrement &) const;

// Other
  boost::shared_ptr<const VariablesQG> variables() const {return varin_;}

  int & toFortran() {return keyOperWind_;}
  const int & toFortran() const {return keyOperWind_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperWind_;
  boost::shared_ptr<const VariablesQG> varin_;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_OBSWINDTLAD_H_
