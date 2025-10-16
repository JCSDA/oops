/*
 * (C) Copyright 2025- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "oops/coupled/GeometryCoupled.h"
#include "oops/coupled/ModelAuxIncrementCoupled.h"
#include "oops/interface/ModelAuxCovariance.h"
#include "oops/util/Printable.h"

namespace oops {

/// Coupled implementation of ModelAuxCovariance
template <typename MODEL1, typename MODEL2>
class ModelAuxCovarianceCoupled : public util::Printable {
  typedef GeometryCoupled<MODEL1, MODEL2>  GeometryCoupled_;
  typedef AuxCoupledModel<MODEL1, MODEL2>  ModelAuxControlCoupled_;
  typedef ModelAuxIncrementCoupled<MODEL1, MODEL2> ModelAuxIncrementCoupled_;
  typedef ModelAuxCovariance<MODEL1>       ModelAuxCovariance1;
  typedef ModelAuxCovariance<MODEL2>       ModelAuxCovariance2;

 public:
  static const std::string classname() {return "oops::ModelAuxCovarianceCoupled";}

  ModelAuxCovarianceCoupled(const eckit::Configuration &, const GeometryCoupled_ &);

  ~ModelAuxCovarianceCoupled() = default;

  void linearize(const ModelAuxControlCoupled_ &, const GeometryCoupled_ &);
  void multiply(const ModelAuxIncrementCoupled_ &, ModelAuxIncrementCoupled_ &) const;
  void inverseMultiply(const ModelAuxIncrementCoupled_ &, ModelAuxIncrementCoupled_ &) const;
  void randomize(ModelAuxIncrementCoupled_ &) const;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<ModelAuxCovariance1> cov1_;
  std::unique_ptr<ModelAuxCovariance2> cov2_;
};

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
ModelAuxCovarianceCoupled<MODEL1, MODEL2>::ModelAuxCovarianceCoupled(
    const eckit::Configuration & conf,
    const GeometryCoupled_ & resol)
  : cov1_(), cov2_() {
  Log::trace() << "ModelAuxCovarianceCoupled::ModelAuxCovarianceCoupled starting" << std::endl;
  util::Timer timer(classname(), "ModelAuxCovarianceCoupled");
  cov1_.reset(new ModelAuxCovariance1(conf, resol.geometry1()));
  cov2_.reset(new ModelAuxCovariance2(conf, resol.geometry2()));
  Log::trace() << "ModelAuxCovarianceCoupled::ModelAuxCovarianceCoupled done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void ModelAuxCovarianceCoupled<MODEL1, MODEL2>::linearize(
    const ModelAuxControlCoupled_ & xx, const GeometryCoupled_ & resol) {
  Log::trace() << "ModelAuxCovarianceCoupled::linearize starting" << std::endl;
  util::Timer timer(classname(), "linearize");
  cov1_->linearize(xx.aux1(), resol.geometry1());
  cov2_->linearize(xx.aux2(), resol.geometry2());
  Log::trace() << "ModelAuxCovarianceCoupled::linearize done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void ModelAuxCovarianceCoupled<MODEL1, MODEL2>::multiply(
    const ModelAuxIncrementCoupled_ & dx1, ModelAuxIncrementCoupled_ & dx2) const {
  Log::trace() << "ModelAuxCovarianceCoupled::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  cov1_->multiply(dx1.modelauxincrement1(), dx2.modelauxincrement1());
  cov2_->multiply(dx1.modelauxincrement2(), dx2.modelauxincrement2());
  Log::trace() << "ModelAuxCovarianceCoupled::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void ModelAuxCovarianceCoupled<MODEL1, MODEL2>::inverseMultiply(
    const ModelAuxIncrementCoupled_ & dx1, ModelAuxIncrementCoupled_ & dx2) const {
  Log::trace() << "ModelAuxCovarianceCoupled::inverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "inverseMultiply");
  cov1_->inverseMultiply(dx1.modelauxincrement1(), dx2.modelauxincrement1());
  cov2_->inverseMultiply(dx1.modelauxincrement2(), dx2.modelauxincrement2());
  Log::trace() << "ModelAuxCovarianceCoupled::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void ModelAuxCovarianceCoupled<MODEL1, MODEL2>::randomize(
    ModelAuxIncrementCoupled_ & dx) const {
  Log::trace() << "ModelAuxCovarianceCoupled::randomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");
  cov1_->randomize(dx.modelauxincrement1());
  cov2_->randomize(dx.modelauxincrement2());
  Log::trace() << "ModelAuxCovarianceCoupled::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void ModelAuxCovarianceCoupled<MODEL1, MODEL2>::print(std::ostream & os) const {
  Log::trace() << "ModelAuxCovarianceCoupled::print starting" << std::endl;
  os << "ModelAuxCovarianceCoupled:" << std::endl;
  if (cov1_) {
    os << *cov1_ << std::endl;
  }
  if (cov2_) {
    os << *cov2_ << std::endl;
  }
  Log::trace() << "ModelAuxCovarianceCoupled::print done" << std::endl;
}

}  // namespace oops
