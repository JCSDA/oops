/*
 * (C) Copyright 2025 UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef LORENZ95_OBSERRORDIAGONAL95_H_
#define LORENZ95_OBSERRORDIAGONAL95_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"

#include "lorenz95/ObsTable.h"
#include "lorenz95/ObsVec1D.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace lorenz95 {

  // Forward declarations
  struct L95ObsTraits;
  class ObsTable;
  class ObsVec1D;

class ObsErrorDiagonal95 : public util::Printable,
                           private boost::noncopyable,
                           private util::ObjectCounter<ObsErrorDiagonal95> {
 public:
  static const std::string classname() {return "lorenz95::ObsErrorDiagonal95";}

  ObsErrorDiagonal95(const eckit::Configuration &, const ObsTable &);
  ~ObsErrorDiagonal95() {}

  void multiply(ObsVec1D &) const;

  void inverseMultiply(ObsVec1D &) const;

  void update(const ObsVec1D &);

  void randomize(ObsVec1D &) const;

  void save(const std::string &) const;

  double getRMSE() const;

  std::unique_ptr<ObsVec1D> getObsErrors() const;

  std::unique_ptr<ObsVec1D> getInverseVariance() const;

 private:
  void print(std::ostream &) const override;

  void randomizeWithoutZeroEnsembleMean(ObsVec1D &) const;

  void randomizeWithZeroEnsembleMean(ObsVec1D &) const;

  ObsVec1D stddev_;
  ObsVec1D inverseVariance_;

  double pert_;
  int member_;
  int numberOfMembers_;
  bool zeroMeanPert_;
};

}  // namespace lorenz95

#endif  // LORENZ95_OBSERRORDIAGONAL95_H_
