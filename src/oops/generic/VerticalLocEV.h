/*
 * (C) Copyright 2020 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_GENERIC_VERTICALLOCEV_H_
#define OOPS_GENERIC_VERTICALLOCEV_H_

#include <Eigen/Dense>

#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/assimilation/Increment4D.h"
#include "oops/base/IncrementEnsemble4D.h"
#include "oops/base/LocalIncrement.h"
#include "oops/generic/gc99.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

namespace oops {

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
};

// ----------------------------------------------------------------------------
template<typename MODEL>
class VerticalLocEV: public util::Printable,
                     private util::ObjectCounter<VerticalLocEV<MODEL>> {
  typedef Geometry<MODEL>            Geometry_;
  typedef GeometryIterator<MODEL>    GeometryIterator_;
  typedef Increment4D<MODEL>         Increment4D_;
  typedef IncrementEnsemble4D<MODEL> IncrementEnsemble4D_;

 public:
  static const std::string classname() {return "oops::VerticalLocEV";}

  VerticalLocEV(const Geometry_ & , const eckit::Configuration &);

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

// replicate eigen vectors for state vars
  std::vector<double> replicateEigenVector(const oops::Variables &,
                                           size_t) const;

 private:
  void print(std::ostream &) const {}
  VerticalLocalizationParameters options_;
  Eigen::MatrixXd Evecs_;
  Eigen::VectorXd Evals_;
  size_t neig_;
};
// -----------------------------------------------------------------------------
template<typename MODEL>
  VerticalLocEV<MODEL>::VerticalLocEV(const Geometry_ & geom,
                               const eckit::Configuration & conf) {
    // read vertical localization configuration
    options_.deserialize(conf);

    // compute vertical corrleation matrix fo nLevs
    Eigen::MatrixXd cov = computeCorrMatrix(geom);
    // compute truncated correlation matrix
    computeCorrMatrixEvec(cov);
    neig_ = truncateEvecs();
  }

// -------------------------------------------------------------------------------------------------
template<typename MODEL>
  Eigen::MatrixXd VerticalLocEV<MODEL>::computeCorrMatrix(const Geometry_ & geom) {
    std::string locUnits = options_.VertLocUnits;
    oops::Log::debug() << "locUnits: " << locUnits << std::endl;
    std::vector<double> vCoord = geom.verticalCoord(locUnits);
    size_t nlevs = vCoord.size();

    // compute vertical correlations and eigen vectors
    Eigen::MatrixXd cov(nlevs, nlevs);
    for (size_t jj=0; jj < nlevs; ++jj) {
      for (size_t ii=jj; ii < nlevs; ++ii) {
        cov(ii, jj) = oops::gc99(std::abs(vCoord[jj]-vCoord[ii])/options_.VertLocDist);
      }
    }
    oops::Log::debug() << "corr: " << cov << std::endl;

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
  // make reference solution
  Eigen::MatrixXd cov = computeCorrMatrix(geom);
  // only the lower traingle of cov is stored
  // transform to full matrix
  Eigen::MatrixXd covOld = cov+cov.transpose();
  covOld.diagonal() = covOld.diagonal()/2;

  // recreate covNew from stored Evecs
  Eigen::MatrixXd covNew = Evecs_*Evecs_.transpose();

  oops::Log::debug() << "Trace covNew " << covNew.diagonal().sum() << std::endl;
  oops::Log::debug() << "Trace cov " << covOld.diagonal().sum() << std::endl;

  oops::Log::debug() << "norm(covNew) " << covNew.norm() << std::endl;
  oops::Log::debug() << "norm(cov) " << covOld.norm() << std::endl;

  covNew = covNew.array() - covOld.array();
  oops::Log::debug() << "norm(covNew-cov) " << covNew.norm() << std::endl;

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
  oops::Variables vars = incr[0].variables();

  for (size_t ieig=0; ieig < neig_; ++ieig) {
    std::vector<double> EvecRepl = replicateEigenVector(vars, ieig);
    // modulate an increment
    for (size_t itime=0; itime < incr.size(); ++itime) {
      const Geometry_ & geom = incr[itime].geometry();
      for (GeometryIterator_ gpi = geom.begin(); gpi != geom.end(); ++gpi) {
        oops::LocalIncrement gp = incr[itime].getLocal(gpi);
        gp *= EvecRepl;
        incrsOut[ieig][itime].setLocal(gp, gpi);
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

  oops::Variables vars = incrs[0][0].variables();
  size_t nv = Evecs_.rows()*vars.size();
  size_t nens = incrs.size();
  Eigen::MatrixXd Z(nv, neig_*nens);
  std::vector<double> etmp2(nv);

  size_t ii = 0;
  for (size_t iens=0; iens < nens; ++iens) {
    oops::LocalIncrement gp = incrs[iens][itime].getLocal(gi);
    etmp2 = gp.getVals();
    for (size_t ieig=0; ieig < neig_; ++ieig) {
      std::vector<double> EvecRepl = replicateEigenVector(vars, ieig);
      // modulate and assign
      for (size_t iv=0; iv < nv; ++iv) {
         Z(iv, ii) = EvecRepl[iv]*etmp2[iv];
      }
      ii += 1;
    }
  }
  return Z;
}
template<typename MODEL>
std::vector<double> VerticalLocEV<MODEL>::replicateEigenVector(
                        const oops::Variables & v, size_t ieig) const {
  // populate eigen vector for all variables nLevs*nVars
  std::vector<double> etmp(Evecs_.rows()*v.size());

  size_t ii = 0;
  for (unsigned ivar=0; ivar < v.size(); ++ivar) {
    for (unsigned ilev=0; ilev < Evecs_.rows(); ++ilev) {
      etmp[ii] = Evecs_(ilev, ieig);
      ii += 1;
    }
  }
  return etmp;
}

}  // namespace oops

#endif  // OOPS_GENERIC_VERTICALLOCEV_H_
