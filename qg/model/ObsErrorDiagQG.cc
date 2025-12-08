/*
 * (C) Copyright 2025 UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/ObsErrorDiagQG.h"

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

void ObsErrorDiagQG::print(std::ostream & os) const {
  os << "Diagonal QG observation error covariance" << std::endl << stddev_;
}

}  // end namespace qg
