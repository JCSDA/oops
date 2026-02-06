/*
 * (C) Copyright 2025 UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <cmath>
#include <sstream>
#include <vector>

#include "lorenz95/ObsErrorDiagonal95.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

namespace lorenz95 {

ObsErrorDiagonal95::ObsErrorDiagonal95(const eckit::Configuration & conf, const ObsTable & ot)
                                    : stddev_(ot, "ObsError"), inverseVariance_(ot, ""),
                                      pert_(conf.getDouble("obs perturbations amplitude", 1.0)),
                                      member_(conf.getInt("member", 1)),
                                      numberOfMembers_(conf.getInt("number of members", 1)),
                                      zeroMeanPert_(conf.getBool("zero-mean perturbations", false))
  {
    oops::Log::trace() << this->classname() << " started" << std::endl;
    inverseVariance_ = stddev_;
    inverseVariance_ *= stddev_;
    inverseVariance_.invert();
    oops::Log::trace() << this->classname() << " constructed" << std::endl;
  }

void ObsErrorDiagonal95::multiply(ObsVec1D & dy) const {
  oops::Log::trace() << "ObsErrorDiagonal95::multiply started" << std::endl;
  dy /= inverseVariance_;
  oops::Log::trace() << "ObsErrorDiagonal95::multiply finished" << std::endl;
}

void ObsErrorDiagonal95::inverseMultiply(ObsVec1D & dy) const {
  oops::Log::trace() << "ObsErrorDiagonal95::inverseMultiply started" << std::endl;
  dy *= inverseVariance_;
  oops::Log::trace() << "ObsErrorDiagonal95::inverseMultiply finished" << std::endl;
}

void ObsErrorDiagonal95::update(const ObsVec1D & obsError) {
  stddev_ = obsError;
  inverseVariance_ = stddev_;
  inverseVariance_ *= stddev_;
  inverseVariance_.invert();
  oops::Log::trace() << this->classname() << " covariance updated " << stddev_.nobs() << std::endl;
}

void ObsErrorDiagonal95::randomize(ObsVec1D & dy) const {
  if (this->zeroMeanPert_)
    randomizeWithZeroEnsembleMean(dy);
  else
    randomizeWithoutZeroEnsembleMean(dy);
}

void ObsErrorDiagonal95::save(const std::string & name) const {
  stddev_.save(name);
}

double ObsErrorDiagonal95::getRMSE() const {
  return stddev_.rms();
}

std::unique_ptr<ObsVec1D> ObsErrorDiagonal95::getObsErrors() const {
  return std::make_unique<ObsVec1D>(stddev_);
}

std::unique_ptr<ObsVec1D> ObsErrorDiagonal95::getInverseVariance() const {
  return std::make_unique<ObsVec1D>(inverseVariance_);
}

void ObsErrorDiagonal95::localize(ObsVec1D & locvector) const {
  oops::Log::trace() << "lorenz95::ObsErrorDiagonal95::localize start" << std::endl;

  const double missing = util::missingValue<double>();

  assert(locvector.size() == inverseVariance_.size());
  std::vector<double> localinvvar;
  for (size_t jj = 0; jj < locvector.size(); ++jj) {
    if (locvector[jj] != missing && locvector[jj] <= 0) {
      throw eckit::BadValue("Localization weights must be positive. Use "
                            "oops::util::missingValue<double>() to indicate "
                            "an observation with a weight of zero.");
    }
    if (locvector[jj] != missing && stddev_[jj] != missing) {
      localinvvar.push_back(locvector[jj] * std::pow(stddev_[jj], -2.0));
    }
  }
  local_inverseVariance_ =
      Eigen::Map<Eigen::VectorXd>(localinvvar.data(), localinvvar.size());
}

Eigen::MatrixXf ObsErrorDiagonal95::localInverseMultiply(const Eigen::MatrixXf & zz) const {
  Eigen::MatrixXf zzRinv(zz.rows(), zz.cols());
  for (int ii = 0; ii < zz.rows(); ++ii) {
    zzRinv(ii, Eigen::all) = zz(ii, Eigen::all)
        .cwiseProduct(local_inverseVariance_.cast<float>().transpose());
  }
  return zzRinv;
}

int ObsErrorDiagonal95::localDim() const {
  return local_inverseVariance_.size();
}

void ObsErrorDiagonal95::print(std::ostream & os) const {
  os << "Diagonal l95 observation error covariance" << std::endl << stddev_;
}

void ObsErrorDiagonal95::randomizeWithoutZeroEnsembleMean(ObsVec1D & dy) const {
  dy.random();
  dy *= stddev_;
  dy *= this->pert_;
}

void ObsErrorDiagonal95::randomizeWithZeroEnsembleMean(ObsVec1D & dy) const {
  ObsVec1D perturbation(dy);
  ObsVec1D sum(dy);
  sum.zero();

  // Generate initial independent perturbations for all ensemble members.
  // Calculate their sum and store this member's perturbations in 'dy'.
  for (int member = 1; member <= this->numberOfMembers_; ++member) {
    perturbation.random();
    sum += perturbation;
    if (member == this->member_)
      dy = perturbation;
  }

  // Subtract the ensemble mean of perturbations from this member's perturbations.
  dy.axpy(-1.0 / this->numberOfMembers_, sum);

  // Scale perturbations to the requested amplitude.
  dy *= stddev_;
  dy *= std::sqrt(this->numberOfMembers_ / (this->numberOfMembers_ - 1.0)) * this->pert_;
}

}  // end namespace lorenz95
