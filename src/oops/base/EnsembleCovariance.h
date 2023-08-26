/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_ENSEMBLECOVARIANCE_H_
#define OOPS_BASE_ENSEMBLECOVARIANCE_H_

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/system/ResourceUsage.h"

#include "oops/assimilation/GMRESR.h"
#include "oops/base/Geometry.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/Increment4D.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/Localization.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/interface/LinearVariableChange.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Timer.h"

namespace oops {

/// Parameters describing generic ensemble covariances.
template <typename MODEL>
class EnsembleCovarianceParameters : public ModelSpaceCovarianceParametersBase<MODEL> {
  OOPS_CONCRETE_PARAMETERS(EnsembleCovarianceParameters,
                           ModelSpaceCovarianceParametersBase<MODEL>)

  typedef typename Increment<MODEL>::ReadParameters_ IncrementReadParameters_;
  typedef typename LinearVariableChange<MODEL>::Parameters_ LinearVarChangeParameters_;
 public:
  /// Parameters for ensemble of increments used in the covariances.
  IncrementEnsembleFromStatesParameters<MODEL> ensemble{this};
  OptionalParameter<IncrementReadParameters_> inflationField{"inflation field",
                   "inflation field (local)", this};
  Parameter<double> inflationValue{"inflation value", "inflation value (global)", 1.0, this};
  OptionalParameter<LinearVarChangeParameters_> ensTrans{"ensemble transform",
                   "ensemble transform: inverse is applied to ensemble members, "
                   "and forward/adjoint around the localized covariance matrix: "
                   "T in the covariance matrix T ( (Tinv X) (Tinv X)t o L ) Tt ", this};
  oops::OptionalParameter<eckit::LocalConfiguration> localization{"localization",
                         "localization applied to ensemble covariances", this};
};

/// Generic ensemble based model space error covariance.

// -----------------------------------------------------------------------------
template <typename MODEL>
class EnsembleCovariance : public ModelSpaceCovarianceBase<MODEL>,
                           private util::ObjectCounter<EnsembleCovariance<MODEL>> {
  typedef Geometry<MODEL>                           Geometry_;
  typedef Increment<MODEL>                          Increment_;
  typedef Increment4D<MODEL>                        Increment4D_;
  typedef LinearVariableChange<MODEL>               LinearVariableChange_;
  typedef Localization<MODEL>                       Localization_;
  typedef State4D<MODEL>                            State4D_;

 public:
  typedef EnsembleCovarianceParameters<MODEL> Parameters_;

  static const std::string classname() {return "oops::EnsembleCovariance";}

  EnsembleCovariance(const Geometry_ &, const Variables &,
                     const Parameters_ &, const State4D_ &, const State4D_ &);
  ~EnsembleCovariance();

 private:
  void doRandomize(Increment4D_ &) const override;
  void doMultiply(const Increment4D_ &, Increment4D_ &) const override;
  void doInverseMultiply(const Increment4D_ &, Increment4D_ &) const override;

  std::unique_ptr<IncrementSet<MODEL>> ens_;
  std::unique_ptr<LinearVariableChange_> ensTrans_;
  Variables ensTransInputVars_;
  Variables ensTransOutputVars_;
  std::unique_ptr<Localization_> loc_;
  int seed_ = 7;  // For reproducibility
};

// -----------------------------------------------------------------------------
template<typename MODEL>
EnsembleCovariance<MODEL>::EnsembleCovariance(const Geometry_ & resol, const Variables & vars,
                                              const Parameters_ & params,
                                              const State4D_ & xb, const State4D_ & fg)
  : ModelSpaceCovarianceBase<MODEL>(resol, params, xb, fg), ens_(),
    ensTransInputVars_(vars), ensTransOutputVars_(vars), loc_()
{
  Log::trace() << "EnsembleCovariance::EnsembleCovariance start" << std::endl;
  util::Timer timer("oops::Covariance", "EnsembleCovariance");
  size_t init = eckit::system::ResourceUsage().maxResidentSetSize();

  eckit::LocalConfiguration conf = params.toConfiguration();

  // Create ensemble
  StateSet<MODEL> tmp(resol, conf, xb.commTime());
  if (tmp.ens_size() < 2) {
    throw eckit::BadParameter("Not enough ensemble members provided for ensemble "
                              "covariances (at least 2 required)", Here());
  }
  ens_.reset(new IncrementSet<MODEL>(resol, vars, tmp, true));
  *ens_ -= ens_->ens_mean();

  // Setup ensemble transform
  if (conf.has("ensemble transform")) {
    // Create ensemble transform
    const auto & ensTransParams = *params.ensTrans.value();
    ensTrans_.reset(new LinearVariableChange_(resol, ensTransParams));

    // Define localization variables as ensemble transform input variables
    // If missing, default is vars (ensemble variables)
    if (ensTransParams.inputVariables.value() != boost::none) {
      ensTransInputVars_ = *ensTransParams.inputVariables.value();
    }

    // If present, check that ensemble transform output variables are
    // a subset of ensemble variables
    // If missing, default is vars (ensemble variables)
    if (ensTransParams.outputVariables.value() != boost::none) {
      ensTransOutputVars_ = *ensTransParams.outputVariables.value();
      ASSERT(ensTransOutputVars_ <= vars);
    }

    // Set trajectory
    ensTrans_->changeVarTraj(fg[0], ensTransOutputVars_);   // change var per time????????
  }

  if (conf.has("inflation value") || conf.has("inflation field")) {
    // Read inflation field
    std::unique_ptr<Increment_> inflationField;
    if (conf.has("inflation field")) {
      inflationField.reset(new Increment_(resol, vars, xb[0].validTime()));
      inflationField->read(*params.inflationField.value());
    }

    // Get inflation value
    const bool inflate = conf.has("inflation value");
    const double inflationValue = conf.getDouble("inflation value", 1.0);

    // Loop over members
    for (size_t jt = 0; jt < ens_->local_time_size(); ++jt) {
      for (size_t jm = 0; jm < ens_->local_ens_size(); ++jm) {
        // Apply local inflation
        if (inflationField) {
          (*ens_)(jt, jm).schur_product_with(*inflationField);
        }

        // Apply global inflation
        if (inflate) (*ens_)(jt, jm) *= inflationValue;

        // Apply ensemble transform inverse
        if (ensTrans_) {
          ensTrans_->changeVarInverseTL((*ens_)(jt, jm), ensTransInputVars_);
        }
      }
    }
  }

  // Create localization
  if (conf.has("localization")) {
    eckit::LocalConfiguration locconf(conf, "localization");
    locconf.set("date", xb[0].validTime().toString());
    locconf.set("time rank", resol.timeComm().rank());
    oops::Variables locVars(vars);
    if (ensTrans_) {
      locVars -= ensTransOutputVars_;
      locVars += ensTransInputVars_;
    }
    loc_.reset(new Localization_(resol, locVars, locconf));
  }

  size_t current = eckit::system::ResourceUsage().maxResidentSetSize();
  this->setObjectSize(current - init);
  Log::trace() << "EnsembleCovariance::EnsembleCovariance done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
EnsembleCovariance<MODEL>::~EnsembleCovariance() {
  Log::trace() << "EnsembleCovariance destructed." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void EnsembleCovariance<MODEL>::doRandomize(Increment4D_ & dx) const {
  ASSERT(ens_->ens_size() == ens_->local_ens_size());
  dx.zero();

  if (loc_) {
    // Localized covariance matrix
    Increment4D_ tmp(dx);
    for (size_t jm = 0; jm < ens_->ens_size(); ++jm) {
      loc_->randomize(tmp);
      for (size_t jt = 0; jt < dx.local_time_size(); ++jt) {
        tmp[jt].schur_product_with((*ens_)(jt, jm));
      }
      dx += tmp;
    }
  } else {
    // Raw covariance matrix
    util::NormalDistribution<double> normalDist(ens_->ens_size(), 0.0, 1.0, seed_);
    for (size_t jm = 0; jm < ens_->ens_size(); ++jm) {
      for (size_t jt = 0; jt < dx.local_time_size(); ++jt) {
        dx[jt].axpy(normalDist[jm], (*ens_)(jt, jm));
      }
    }
  }

  // Ensemble transform
  if (ensTrans_) {
    for (size_t jt = 0; jt < dx.local_time_size(); ++jt) {
      ensTrans_->changeVarTL(dx[jt], ensTransOutputVars_);
    }
  }

  // Normalization
  dx *= 1.0/sqrt(static_cast<double>(ens_->size()-1));
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void EnsembleCovariance<MODEL>::doMultiply(const Increment4D_ & dxi, Increment4D_ & dxo) const {
  ASSERT(ens_->ens_size() == ens_->local_ens_size());
  Increment4D_ dxiTmp(dxi);

  // Ensemble transform adjoint
  if (ensTrans_) {
    for (size_t jt = 0; jt < ens_->local_time_size(); ++jt) {
      ensTrans_->changeVarAD(dxiTmp[jt], ensTransInputVars_);
    }
  }

  Increment4D_ dxoTmp(dxiTmp);
  dxoTmp.zero();

  // Localization
  if (loc_) {
    // Localized covariance matrix
    for (size_t jm = 0; jm < ens_->ens_size(); ++jm) {
      Increment4D_ dx(dxiTmp);
      for (size_t jt = 0; jt < ens_->local_time_size(); ++jt) {
        dx[jt].schur_product_with((*ens_)(jt, jm));
      }
      loc_->multiply(dx);
      for (size_t jt = 0; jt < ens_->local_time_size(); ++jt) {
        dx[jt].schur_product_with((*ens_)(jt, jm));
      }
      dxoTmp += dx;
    }
  } else {
    for (size_t jm = 0; jm < ens_->ens_size(); ++jm) {
      double wgt = 0;
      for (size_t jt = 0; jt < ens_->local_time_size(); ++jt) {
        wgt += dxiTmp[jt].dot_product_with((*ens_)(jt, jm));
      }
      dxi.geometry().timeComm().allReduceInPlace(wgt, eckit::mpi::sum());
      for (size_t jt = 0; jt < ens_->local_time_size(); ++jt) {
        dxoTmp[jt].axpy(wgt, (*ens_)(jt, jm), false);
      }
    }
  }

  // Ensemble transform
  if (ensTrans_) {
    for (size_t jt = 0; jt < ens_->local_time_size(); ++jt) {
      ensTrans_->changeVarTL(dxoTmp[jt], ensTransOutputVars_);
    }
  }

  // Normalization
  dxoTmp *= 1.0/static_cast<double>(ens_->ens_size()-1);

  // Copy to output Increment_
  dxo = dxoTmp;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void EnsembleCovariance<MODEL>::doInverseMultiply(const Increment4D_ & dxi,
                                                  Increment4D_ & dxo) const {
  IdentityMatrix<Increment4D_> Id;
  dxo.zero();
  GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_ENSEMBLECOVARIANCE_H_
