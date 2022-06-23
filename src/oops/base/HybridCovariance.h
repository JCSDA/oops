/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_HYBRIDCOVARIANCE_H_
#define OOPS_BASE_HYBRIDCOVARIANCE_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "oops/assimilation/GMRESR.h"
#include "oops/base/Geometry.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/Increment.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace oops {

/// Generic hybrid static-ensemble model space error covariance.

// -----------------------------------------------------------------------------
template <typename MODEL>
class HybridCovariance : public ModelSpaceCovarianceBase<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
  HybridCovariance(const Geometry_ &, const Variables &,
                   const eckit::Configuration &, const State_ &, const State_ &);
  ~HybridCovariance();

 private:
  void doRandomize(Increment_ &) const override;
  void doMultiply(const Increment_ &, Increment_ &) const override;
  void doInverseMultiply(const Increment_ &, Increment_ &) const override;

  std::vector< std::unique_ptr< ModelSpaceCovarianceBase<MODEL> > > Bcomponents_;
  std::vector< std::string > weightTypes_;
  std::vector< double > valueWeights_;
  std::vector< Increment_ > incrementWeightsSqrt_;
};

// =============================================================================

/// Constructor, destructor
// -----------------------------------------------------------------------------
template<typename MODEL>
HybridCovariance<MODEL>::HybridCovariance(const Geometry_ & resol, const Variables & vars,
                                          const eckit::Configuration & config,
                                          const State_ & xb, const State_ & fg)
  : ModelSpaceCovarianceBase<MODEL>(resol, config, xb, fg)
{
  Log::trace() << "HybridCovariance::HybridCovariance start" << std::endl;
  util::Timer timer("oops::Covariance", "HybridCovariance");
  std::vector<eckit::LocalConfiguration> confs;
  config.get("components", confs);
  for (const auto & conf : confs) {
    // B component
    const eckit::LocalConfiguration covarConf(conf, "covariance");
    std::unique_ptr< ModelSpaceCovarianceBase<MODEL> > B(
        CovarianceFactory<MODEL>::create(resol, vars, covarConf, xb, fg));
    Bcomponents_.push_back(std::move(B));

    // Weight
    const eckit::LocalConfiguration weightConf(conf, "weight");
    if (weightConf.has("value")) {
      // Scalar weight provided in the configuration
      weightTypes_.push_back("value");
      const double valueWeight = weightConf.getDouble("value");
      ASSERT(valueWeight >= 0.0);
      valueWeights_.push_back(valueWeight);
    } else {
      // 3D weight read from a file
      weightTypes_.push_back("increment");
      Increment_ weight(resol, vars, xb.validTime());
      weight.read(weightConf);

      // Compute weight square-root
      util::FieldSetSqrt(weight.fieldSet());
      weight.synchronizeFields();
      incrementWeightsSqrt_.push_back(weight);
    }
  }
  Log::trace() << "HybridCovariance::HybridCovariance done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
HybridCovariance<MODEL>::~HybridCovariance() {
  Log::trace() << "HybridCovariance destructed" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void HybridCovariance<MODEL>::doMultiply(const Increment_ & dxi, Increment_ & dxo) const {
  dxo.zero();
  Increment_ tmp(dxo);
  int valueIndex = 0;
  int incrementIndex = 0;
  for (size_t jcomp = 0; jcomp < Bcomponents_.size(); ++jcomp) {
     if (weightTypes_[jcomp] == "increment") {
        tmp.schur_product_with(incrementWeightsSqrt_[incrementIndex]);
     }
     Bcomponents_[jcomp]->multiply(dxi, tmp);
     if (weightTypes_[jcomp] == "value") {
        tmp *= valueWeights_[valueIndex];
        valueIndex += 1;
     }
     if (weightTypes_[jcomp] == "increment") {
        tmp.schur_product_with(incrementWeightsSqrt_[incrementIndex]);
        incrementIndex += 1;
     }
     dxo += tmp;
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void HybridCovariance<MODEL>::doInverseMultiply(const Increment_ & dxi, Increment_ & dxo) const {
  IdentityMatrix<Increment_> Id;
  dxo.zero();
  GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void HybridCovariance<MODEL>::doRandomize(Increment_ & dx) const {
  dx.zero();
  Increment_ tmp(dx);
  int valueIndex = 0;
  int incrementIndex = 0;
  for (size_t jcomp = 0; jcomp < Bcomponents_.size(); ++jcomp) {
     Bcomponents_[jcomp]->randomize(tmp);
     if (weightTypes_[jcomp] == "value") {
        tmp *= std::sqrt(valueWeights_[valueIndex]);
        valueIndex += 1;
     }
     if (weightTypes_[jcomp] == "increment") {
        tmp.schur_product_with(incrementWeightsSqrt_[incrementIndex]);
        incrementIndex += 1;
     }
     dx += tmp;
  }
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_HYBRIDCOVARIANCE_H_
