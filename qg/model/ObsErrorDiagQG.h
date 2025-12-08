/*
 * (C) Copyright 2025 UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef QG_MODEL_OBSERRORDIAGQG_H_
#define QG_MODEL_OBSERRORDIAGQG_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"

#include "model/ObsSpaceQG.h"
#include "model/ObsVecQG.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace qg {
class ObsSpaceQG;
class ObsVecQG;
struct QgObsTraits;

class ObsErrorDiagQG : public util::Printable,
                       private boost::noncopyable,
                       private util::ObjectCounter<ObsErrorDiagQG> {
 public:
  static const std::string classname() {return "qg::ObsErrorDiagQG";}

  ObsErrorDiagQG(const eckit::Configuration &, const ObsSpaceQG &);
  ~ObsErrorDiagQG() {}

  void multiply(ObsVecQG &) const;

  void inverseMultiply(ObsVecQG &) const;

  void update(const ObsVecQG &);

  void randomize(ObsVecQG &) const;

  void save(const std::string &) const;

  double getRMSE() const;

  std::unique_ptr<ObsVecQG> getObsErrors() const;

  std::unique_ptr<ObsVecQG> getInverseVariance() const;

 private:
  void print(std::ostream &) const override;

  ObsVecQG stddev_;
  ObsVecQG inverseVariance_;
  double pert_;
};

}  // end namespace qg

#endif  // QG_MODEL_OBSERRORDIAGQG_H_
