/*
 * (C) Copyright 2020-2021 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_GENERIC_VERTICALLOCEV_H_
#define OOPS_GENERIC_VERTICALLOCEV_H_

#include <Eigen/Dense>

#include <memory>
#include <string>
#include <vector>

#include <boost/make_unique.hpp>

#include "eckit/config/Configuration.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment4D.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/IncrementEnsemble4D.h"
#include "oops/base/LocalIncrement.h"
#include "oops/generic/gc99.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

namespace oops {
  class Variables;

/// Parameters for vertical localization
/*!
 * See Lei 2018 JAMES for more details
 *
 * Lei, L., Whitaker, J. S., & Bishop, C. ( 2018). Improving assimilation
 * of radiance observations by implementing model space localization in an
 * ensemble Kalman filter. Journal of Advances in Modeling Earth Systems, 10,
 * 3221– 3232. https://doi.org/10.1029/2018MS001468
 */
class VerticalLocalizationParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(VerticalLocalizationParameters, Parameters)
 public:
  // fraction of the variance retained after the eigen spectrum
  // of the vertical localization function is truncated
  // 1 -- retain all eigen vectors
  // 0 -- retain the first eigen vector
  Parameter<double> VertLocToll{"fraction of retained variance", 1.0, this};
  // localization distance at which Gaspari-Cohn = 0
  RequiredParameter<double> VertLocDist{"lengthscale", this};
  // localization distance at which Gaspari-Cohn = 0
  RequiredParameter<std::string> VertLocUnits{"lengthscale units", this};
  // write eigen vectors to disk
  Parameter<bool> writeEVs{"write eigen vectors", false, this};
  // read eigen vectors from disk
  Parameter<bool> readEVs{"read eigen vectors", false, this};
};

// ----------------------------------------------------------------------------
template<typename MODEL>
class VerticalLocEV: public util::Printable,
                     private util::ObjectCounter<VerticalLocEV<MODEL>> {
  typedef Geometry<MODEL>            Geometry_;
  typedef GeometryIterator<MODEL>    GeometryIterator_;
  typedef Increment4D<MODEL>         Increment4D_;
  typedef IncrementEnsemble<MODEL>   IncrementEnsemble_;
  typedef IncrementEnsemble4D<MODEL> IncrementEnsemble4D_;
  typedef State<MODEL>               State_;

 public:
  static const std::string classname() {return "oops::VerticalLocEV";}

  VerticalLocEV(const eckit::Configuration &, const State_ &, const Variables &);

// modulate an increment
  void modulateIncrement(const Increment4D_ &, IncrementEnsemble4D_ &) const;

// modulate an incrementEnsemble at a {gridPoint, timeSlice}
  Eigen::MatrixXd modulateIncrement(const IncrementEnsemble4D_ &,
                                    const GeometryIterator_ &, size_t) const;

// returns number of retained eigen modes
  size_t neig() const {return neig_;}

// tests if the truncation and rescaling are done properly
  bool testTruncateEvecs(const Geometry_ &);

 private:
// compute vertical correlations and eigen vectors
  Eigen::MatrixXd computeCorrMatrix(const Geometry_ &);
  void computeCorrMatrixEvec(const Eigen::MatrixXd &);

// truncate Evec sequence
  size_t truncateEvecs();

// replicate eigen vectors for increment vars
  std::vector<double> replicateEigenVector(size_t) const;

// populate 3D increment array of eigen vectors
  void populateIncrementEnsembleWithEVs();
// IO for eigen vectors
  void writeEVsToDisk(const eckit::Configuration &) const;
  void readEVsFromDisk(const Geometry_ &, const util::DateTime &,
                       const eckit::Configuration &);

 private:
  void print(std::ostream &) const {}
  VerticalLocalizationParameters options_;
  Eigen::MatrixXd Evecs_;
  Eigen::VectorXd Evals_;
  size_t neig_;
  std::unique_ptr<IncrementEnsemble_> sqrtVertLoc_;
  // if true store EVs as explicit 3D fields
  // TODO(frolovsa) depriciate 1D storage in the future
  // once we see no issues with always using 3D
  bool EVsStoredAs3D_ = true;  // if true store EVs as explicit 3D fields
  // increment variables for which vertical localization will be computed
  const Variables incvars_;
};
// -----------------------------------------------------------------------------
template<typename MODEL>
  VerticalLocEV<MODEL>::VerticalLocEV(const eckit::Configuration & conf,
                               const State_ & x, const Variables & incvars):
  sqrtVertLoc_(), incvars_(incvars) {
    // read vertical localization configuration
    options_.deserialize(conf);

    if (options_.readEVs) {
      oops::Log::info() << "Reading precomputed vertical localization EVs from disk" << std::endl;
      readEVsFromDisk(x.geometry(), x.validTime(), conf);
      neig_ = sqrtVertLoc_->size();
    } else {
      oops::Log::info() << "Computing vertical localization EVs from scratch" << std::endl;
      // compute vertical corrleation matrix fo nLevs
      Eigen::MatrixXd cov = computeCorrMatrix(x.geometry());
      // compute truncated correlation matrix
      computeCorrMatrixEvec(cov);
      neig_ = truncateEvecs();
      // convert EVs to 3D if needed
      if (EVsStoredAs3D_) {
        sqrtVertLoc_ = boost::make_unique<IncrementEnsemble_>
                        (x.geometry(), incvars_, x.validTime(), neig_);
        populateIncrementEnsembleWithEVs();
      }
    }

    if (options_.writeEVs) { writeEVsToDisk(conf); }
  }

// -----------------------------------------------------------------------------
template<typename MODEL>
  void VerticalLocEV<MODEL>::populateIncrementEnsembleWithEVs() {
    const Geometry_ & geom = (*sqrtVertLoc_)[0].geometry();
    for (size_t ieig=0; ieig < neig_; ++ieig) {
      std::vector<double> EvecRepl = replicateEigenVector(ieig);
      (*sqrtVertLoc_)[ieig].ones();
      for (GeometryIterator_ gpi = geom.begin(); gpi != geom.end(); ++gpi) {
        oops::LocalIncrement gp = (*sqrtVertLoc_)[ieig].getLocal(gpi);
        gp *= EvecRepl;
        (*sqrtVertLoc_)[ieig].setLocal(gp, gpi);
      }
    }
  }

// -----------------------------------------------------------------------------
template<typename MODEL>
  void VerticalLocEV<MODEL>::writeEVsToDisk(const eckit::Configuration & config) const {
    eckit::LocalConfiguration outConfig(config, "list of eigen vectors to write");
    sqrtVertLoc_->write(outConfig);
  }

// -----------------------------------------------------------------------------
template<typename MODEL>
  void VerticalLocEV<MODEL>::readEVsFromDisk(const Geometry_ & geom,
                                             const util::DateTime & tslot,
                                             const eckit::Configuration & config) {
    std::vector<eckit::LocalConfiguration> memberConfig;
    config.get("list of eigen vectors to read", memberConfig);
    sqrtVertLoc_ = boost::make_unique<IncrementEnsemble_>(
                     geom, incvars_, tslot, memberConfig.size());

    // Loop over all ensemble members
    for (size_t jj = 0; jj < memberConfig.size(); ++jj) {
      (*sqrtVertLoc_)[jj].read(memberConfig[jj]);
    }
    // TODO(frolovsa) check that correlation matrix has 1 on the diagnoal.
    // If not, assume that fewer vectors were read and renormalize
  }

// -------------------------------------------------------------------------------------------------
template<typename MODEL>
  Eigen::MatrixXd VerticalLocEV<MODEL>::computeCorrMatrix(const Geometry_ & geom) {
    std::string locUnits = options_.VertLocUnits;
    std::vector<double> vCoord = geom.verticalCoord(locUnits);
    size_t nlevs = vCoord.size();

    // compute vertical correlations and eigen vectors
    Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(nlevs, nlevs);
    for (size_t jj=0; jj < nlevs; ++jj) {
      for (size_t ii=jj; ii < nlevs; ++ii) {
        cov(ii, jj) = oops::gc99(std::abs(vCoord[jj]-vCoord[ii])/options_.VertLocDist);
      }
    }

    return cov;
  }


// -------------------------------------------------------------------------------------------------
template<typename MODEL>
  void VerticalLocEV<MODEL>::computeCorrMatrixEvec(const Eigen::MatrixXd & cov) {
    // compute eigen decomposition of vertical correlation matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(cov);
    Evals_ = es.eigenvalues().real();
    Evecs_ = es.eigenvectors().real();

    // reverse order, largest eigen value first.
    Evals_.reverseInPlace();
    Evecs_.rowwise().reverseInPlace();
  }

// -------------------------------------------------------------------------------------------------
template<typename MODEL>
  size_t VerticalLocEV<MODEL>::truncateEvecs() {
  // truncate Evec sequence to retain "toll" percetage of total variance
    double evalsum = Evals_.sum();

    // assert condition on the trace of the correlation matrix
    if ((evalsum-Evals_.size()) > 1e-5) {
      Log::error() << "VerticalLocEV trace(cov)~=cov.size: trace(cov)=" <<
                       evalsum << "cov.size=" << Evals_.size() << std::endl;
      throw eckit::BadValue("VerticalLocEV trace(cov)~=cov.size");
    }

    // compute number of evals bellow tolerence
    double frac = 0.0;
    size_t neig = 0;
    for (int ii=0; ii < Evals_.size(); ++ii) {
      frac += Evals_[ii];
      neig += 1;
      if (frac/evalsum >= options_.VertLocToll) { break; }
    }
    frac = frac/evalsum;

    oops::Log::info() << "Retaining frac: " << frac << std::endl;
    oops::Log::info() << "Retaining neig vectors in vert. loc.: " << neig << std::endl;

    // compute Eivec(:,i)=Eivec(:,i)*sqrt(Evals(i)/frac)
    // i.e. make Evecs_=sgrt(cov)
    Evals_.array() *= 1/frac;
    for (int ii=0; ii < Evals_.size(); ++ii) {
      Evecs_.array().col(ii) *= std::sqrt(Evals_(ii));
    }

    // truncate to keep neig
    Evecs_.conservativeResize(Eigen::NoChange, neig);
    Evals_.conservativeResize(neig);

    // assert that the trace of the new and old cov is the same
    if ((evalsum-Evals_.sum()) > 1e-5) {
      Log::error() << "VerticalLocEV::truncateEvecs trace(cov)~=trace(covTrunc): " <<
          " trace(cov)=" << evalsum << " trace(covTrunc)=" << Evals_.sum() <<
          " diff=" << (evalsum-Evals_.sum()) << std::endl;
      throw eckit::BadValue("VerticalLocEV::truncateEvecs trace(cov)~=trace(covTrunc)");
    }

    // return number of truncated eigen values
    return neig;
  }


// -----------------------------------------------------------------------------
template<typename MODEL>
bool VerticalLocEV<MODEL>::testTruncateEvecs(const Geometry_ & geom) {
  // only do this test if we are computing EVs from scratch
  // if reading from disk return true
  if (options_.readEVs) {return true;}

  // make reference solution
  Eigen::MatrixXd cov = computeCorrMatrix(geom);
  // only the lower traingle of cov is stored
  // transform to full matrix
  Eigen::MatrixXd covOld = cov+cov.transpose();
  covOld.diagonal() = covOld.diagonal()/2;

  // recreate covNew from stored Evecs
  Eigen::MatrixXd covNew = Evecs_*Evecs_.transpose();

  covNew = covNew.array() - covOld.array();

  bool rv = (covNew.diagonal().sum() - cov.diagonal().sum()) < 1e-5;
  return rv;
}

// -----------------------------------------------------------------------------
template<typename MODEL>
void VerticalLocEV<MODEL>::modulateIncrement(const Increment4D_ & incr,
                                      IncrementEnsemble4D_ & incrsOut) const {
  // modulate an increment incr using Eivec_
  // returns incrsOut

  // create temporary vector that would be used to store the full
  // eigen vector of the correlation matrix
  // at this point Evecs_ stores eigen vector.size()=nLevs
  // we need nLevs*nVars
  // TODO(Issue #812) catch a special case where some variables in the gp are not 3D
  // consider doing it once instead of replicating this for each gp (Issue #813)

  for (size_t ieig=0; ieig < neig_; ++ieig) {
    // modulate an increment
    for (size_t itime=0; itime < incr.size(); ++itime) {
      if (EVsStoredAs3D_) {  // if EVs stored as 3D use oops operation for schur_product
        incrsOut[ieig][itime] = (*sqrtVertLoc_)[ieig];
        incrsOut[ieig][itime].schur_product_with(incr[itime]);
      } else {  // if EV is only available as a 1D column use grid eterator
        std::vector<double> EvecRepl = replicateEigenVector(ieig);
        const Geometry_ & geom = incr[itime].geometry();
        for (GeometryIterator_ gpi = geom.begin(); gpi != geom.end(); ++gpi) {
          oops::LocalIncrement gp = incr[itime].getLocal(gpi);
          gp *= EvecRepl;
          incrsOut[ieig][itime].setLocal(gp, gpi);
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
template<typename MODEL>
Eigen::MatrixXd VerticalLocEV<MODEL>::modulateIncrement(
                                  const IncrementEnsemble4D_ & incrs,
                                  const GeometryIterator_ & gi,
                                  size_t itime) const {
  // modulate an increment at grid point

  size_t nv = 0;
  std::vector<double> EvecRepl;
  if (EVsStoredAs3D_) {
    nv = (*sqrtVertLoc_)[0].getLocal(gi).getVals().size();
  } else {
    EvecRepl = replicateEigenVector(0);
    nv = EvecRepl.size();
  }
  size_t nens = incrs.size();
  Eigen::MatrixXd Z(nv, neig_*nens);
  std::vector<double> etmp2(nv);

  size_t ii = 0;
  for (size_t iens=0; iens < nens; ++iens) {
    etmp2 =  incrs[iens][itime].getLocal(gi).getVals();
    for (size_t ieig=0; ieig < neig_; ++ieig) {
      if (EVsStoredAs3D_) {
        EvecRepl = (*sqrtVertLoc_)[ieig].getLocal(gi).getVals();
      } else {
        EvecRepl = replicateEigenVector(ieig);
      }
      // modulate and assign
      for (size_t iv=0; iv < nv; ++iv) {
         Z(iv, ii) = EvecRepl[iv]*etmp2[iv];
      }
      ii += 1;
    }
  }
  return Z;
}

// -----------------------------------------------------------------------------
template<typename MODEL>
std::vector<double> VerticalLocEV<MODEL>::replicateEigenVector(size_t ieig) const {
  // populate eigen vector for all variables nLevs*nVars
  std::vector<double> etmp(Evecs_.rows()*incvars_.size());

  size_t ii = 0;
  for (unsigned ivar=0; ivar < incvars_.size(); ++ivar) {
    for (unsigned ilev=0; ilev < Evecs_.rows(); ++ilev) {
      etmp[ii] = Evecs_(ilev, ieig);
      ii += 1;
    }
  }
  return etmp;
}

}  // namespace oops

#endif  // OOPS_GENERIC_VERTICALLOCEV_H_
