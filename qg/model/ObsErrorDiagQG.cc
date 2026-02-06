/*
 * (C) Copyright 2025 UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "model/ObsErrorDiagQG.h"

#include "oops/util/missingValues.h"

namespace qg {

ObsErrorDiagQG::ObsErrorDiagQG(const eckit::Configuration & conf, const ObsSpaceQG & ot)
                                        : stddev_(ot, "ObsError"), inverseVariance_(ot, ""),
                                          pert_(conf.getDouble("obs perturbations amplitude", 1.0))
  {
    inverseVariance_ = stddev_;
    inverseVariance_ *= stddev_;
    inverseVariance_.invert();
    oops::Log::trace() << this->classname() << " constructed" << std::endl;
  }

void ObsErrorDiagQG::multiply(ObsVecQG & dy) const {
  dy /= inverseVariance_;
}

void ObsErrorDiagQG::inverseMultiply(ObsVecQG & dy) const {
  dy *= inverseVariance_;
}

void ObsErrorDiagQG::update(const ObsVecQG & obsError) {
  stddev_ = obsError;
  inverseVariance_ = stddev_;
  inverseVariance_ *= stddev_;
  inverseVariance_.invert();
  oops::Log::trace() << this->classname() << " covariance updated " << stddev_.nobs() << std::endl;
}

void ObsErrorDiagQG::randomize(ObsVecQG & dy) const {
  dy.random();
  dy *= stddev_;
  dy *= this->pert_;
}

void ObsErrorDiagQG::save(const std::string & name) const {
  stddev_.save(name);
}

double ObsErrorDiagQG::getRMSE() const {
  return stddev_.rms();
}

std::unique_ptr<ObsVecQG> ObsErrorDiagQG::getObsErrors() const {
  return std::make_unique<ObsVecQG>(stddev_);
}

std::unique_ptr<ObsVecQG> ObsErrorDiagQG::getInverseVariance() const {
  return std::make_unique<ObsVecQG>(inverseVariance_);
}

void ObsErrorDiagQG::localize(ObsVecQG & locvector) const {
  oops::Log::trace() << "qg::ObsErrorDiagQG::localize start" << std::endl;

  const double missing = util::missingValue<double>();

  assert(locvector.size() == inverseVariance_.size());
  std::vector<double> localinvvar;
  std::vector<double> locvectorstd, stddev_std;
  locvector.serialize(locvectorstd);
  stddev_.serialize(stddev_std);
  for (size_t jj = 0; jj < locvector.size(); ++jj) {
    if (locvectorstd[jj] != missing && locvectorstd[jj] <= 0) {
      throw eckit::BadValue("Localization weights must be positive. Use "
                            "oops::util::missingValue<double>() to indicate "
                            "an observation with a weight of zero.");
    }
    if (locvectorstd[jj] != missing && stddev_std[jj] != missing) {
      localinvvar.push_back(locvectorstd[jj] * std::pow(stddev_std[jj], -2.0));
    }
  }
  local_inverseVariance_ =
      Eigen::Map<Eigen::VectorXd>(localinvvar.data(), localinvvar.size());
}

Eigen::MatrixXf ObsErrorDiagQG::localInverseMultiply(const Eigen::MatrixXf & zz) const {
  oops::Log::trace() << "qg::ObsErrorDiagQG::localInverseMultiply start" << std::endl;
  Eigen::MatrixXf zzRinv(zz.rows(), zz.cols());
  for (int ii = 0; ii < zz.rows(); ++ii) {
    zzRinv(ii, Eigen::all) = zz(ii, Eigen::all)
        .cwiseProduct(local_inverseVariance_.cast<float>().transpose());
  }
  return zzRinv;
}

int ObsErrorDiagQG::localDim() const {
  return local_inverseVariance_.size();
}

void ObsErrorDiagQG::print(std::ostream & os) const {
  os << "Diagonal QG observation error covariance" << std::endl << stddev_;
}

}  // end namespace qg
