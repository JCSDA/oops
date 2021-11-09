/*
 * (C) Copyright 2009-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#ifndef OOPS_RUNS_SQRTOFVERTLOC_H_
#define OOPS_RUNS_SQRTOFVERTLOC_H_


#include <algorithm>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/interface/State.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {
/* -----------------------------------------------------------------------------
*  @brief this program computes sqrt of the vertical correlation in B
*  this is done for each grid point by computing the eigen value problem
*  for the correlation matrix computed from random darws from the B matrix
*  the sqrt(B) is output as a truncated sequence of eigen vectors
*  -----------------------------------------------------------------------------
*/
template <typename MODEL> class SqrtOfVertLoc : public Application {
  typedef Geometry<MODEL>            Geometry_;
  typedef GeometryIterator<MODEL>    GeometryIterator_;
  typedef Increment<MODEL>           Increment_;
  typedef IncrementEnsemble<MODEL>   IncrementEnsemble_;
  typedef ModelSpaceCovarianceBase<MODEL> ModelSpaceCovariance_;
  typedef State<MODEL>               State_;

 public:
// -----------------------------------------------------------------------------
  explicit SqrtOfVertLoc(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateCovarFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~SqrtOfVertLoc() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
//  truncation tolerence
    const double TruncationTollerance = fullConfig.getDouble("truncation tolerence", 1.0);

//  Setup geometry
    const eckit::LocalConfiguration geometryConfig(fullConfig, "geometry");
    const Geometry_ geometry(geometryConfig, this->getComm());

//  Setup background state
    const eckit::LocalConfiguration initialConfig(fullConfig, "background");
    const State_ xx(geometry, initialConfig);
    Log::test() << "Background: " << xx << std::endl;

//  Setup variables
    const Variables vars(fullConfig, "perturbed variables");

//  Setup B matrix
    const eckit::LocalConfiguration covar(fullConfig, "background error");
    std::unique_ptr< ModelSpaceCovarianceBase<MODEL> >
      Bmat(CovarianceFactory<MODEL>::create(covar, geometry, vars, xx, xx));

//  Retreave vertical eigen vectors from B
    const size_t samples = fullConfig.getInt("number of random samples");
    IncrementEnsemble_ perts(geometry, vars, xx.validTime(), samples);
    size_t maxNeigOutput = fullConfig.getInt("max neig output", samples);
    size_t truncatedNeig = getVerticalEigenVectors(*Bmat, geometry, perts,
                           TruncationTollerance, maxNeigOutput);

//  Inflate truncated columns of sqrt(B) to account for the truncated spectrum
    // (1) compute trace of the truncated correlation matrix
    Increment_ tmpIncr1(geometry, vars, xx.validTime());
    Increment_ sumOfSquares(geometry, vars, xx.validTime());
    sumOfSquares.zero();
    for (int jm = 0; jm < truncatedNeig; ++jm) {
      tmpIncr1 = perts[jm];
      tmpIncr1.schur_product_with(tmpIncr1);
      sumOfSquares += tmpIncr1;
    }
    oops::Log::test() << "Variance for truncated correlation matrix" << sumOfSquares << std::endl;
    // (2) compute inflation as 1./sqrt(sumOfSquares)
    for (GeometryIterator_ i = geometry.begin(); i != geometry.end(); ++i) {
      oops::LocalIncrement li = sumOfSquares.getLocal(i);
      std::vector<double> doubleVector = li.getVals();
      for (size_t ii = 0; ii < doubleVector.size(); ++ii) {
        if (doubleVector[ii] != 0) { doubleVector[ii] = 1/std::sqrt(doubleVector[ii]); }
      }
      li.setVals(doubleVector);
      sumOfSquares.setLocal(li, i);
    }

//  Output columns of sqrt(B)
    for (int jm = 0; jm < truncatedNeig; ++jm) {
      eckit::LocalConfiguration outConfig(fullConfig, "output");
      outConfig.set("member", jm+1);
      perts[jm].schur_product_with(sumOfSquares);  //  Scale eigen vectors
      perts[jm].write(outConfig);
      Log::test() << "Columns of sqrt(B) " << jm << perts[jm] << std::endl;
    }
    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::SqrtOfVertLoc<" + MODEL::name() + ">";
  }

  /*  @brief samples covariance matrix and computes vertical eigen vectors of B
            perts.size() -- defines size of the randomization
    @param cov   -- covariance matrix
    @param perts -- preallocated structure for IncrementEnsemble
                    on entrence can be anything. just perts.size() is used
                    on exit contains vertical eigen vectors in decreasingorder
    @param TruncationTollerance -- float [0 1] defines trunctation
    @param maxNeigOutput -- maxNeigOutput for output
  */
  size_t getVerticalEigenVectors(const ModelSpaceCovariance_ & cov,
                                 const Geometry_ & geom,
                                 IncrementEnsemble_ & perts,
                                 const double TruncationTollerance,
                                 const size_t maxNeigOutput) const {
//  Mpi communicator
    const eckit::mpi::Comm & mpiComm = geom.getComm();

//  Generate random sample of B
    size_t samples = perts.size();
    for (size_t jm = 0; jm < samples; ++jm) {
      cov.randomize(perts[jm]);
    }

//  Create temp. eigen matrices
    oops::LocalIncrement liTmp = perts[0].getLocal(geom.begin());
    size_t nv = liTmp.getVals().size();  // number of variables in a vertical column
    Eigen::MatrixXd Z(nv, samples);
    Eigen::VectorXd averageEigenSpectrum = Eigen::VectorXd::Zero(std::min(samples, nv));
    Eigen::VectorXd varianceExplained = Eigen::VectorXd::Zero(std::min(samples, nv));
    std::vector<double> dVectorTmp(nv);

//  Loop over the grid points and compute the eigen value problems
    int numberOfPointsOnThisPE = 0;
    for (GeometryIterator_ i = geom.begin(); i != geom.end(); ++i) {
      // make sure we are working with a valid local perturbation
      // e.g. local perts are 0 for ocean points on land
      std::vector<double> doubleVector = perts[0].getLocal(i).getVals();
      double l = std::inner_product(doubleVector.begin(), doubleVector.end(),
                                    doubleVector.begin(), 0.0);
      if ( l == 0 ) {
        // zero-length perturbation, skip eigen value computation
        // do not increment numberOfPointsOnThisPE
        perts.setEigen(Eigen::MatrixXd::Zero(nv, samples), i);
      } else {
        // non-zero perturbation, proceed with eigen value computation

        // populate the eigen matrix
        perts.packEigen(Z, i);

        // compute correlation matrix
        Eigen::MatrixXd Zcentered = Z.colwise() - Z.rowwise().mean();
        Eigen::MatrixXd corr = (Zcentered * Zcentered.adjoint()) /
                                static_cast<double>(Z.cols() - 1);
        Eigen::VectorXd invZstd = corr.diagonal().cwiseSqrt().cwiseInverse();
        invZstd = (invZstd.array().isFinite()).select(invZstd, 0);
        corr = invZstd.asDiagonal()*corr*invZstd.asDiagonal();

        // compute the eigen decomposition of the correlation matrix
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(corr);
        Eigen::MatrixXd Evals = es.eigenvalues().real();
        Eigen::MatrixXd Evecs = es.eigenvectors().real();

        // reverse order, largest eigen value first.
        Evals.reverseInPlace();
        Evecs.rowwise().reverseInPlace();

        // convert to square root
        // compute Evec(:,i)=Evec(:,i)*sqrt(Evals(i))
        // i.e. make Evecs=sgrt(cov)
        for (int ii=0; ii < Evals.size(); ++ii) {
          Evecs.array().col(ii) *= std::sqrt(Evals(ii));
        }

        // assign the result from the eigen object to IncrementEnsemble
        Z.setZero();
        Z.leftCols(std::min(Evecs.cols(), Z.cols())) =
                   Evecs.leftCols(std::min(Evecs.cols(), Z.cols()));
        perts.setEigen(Z, i);

        // increment the eigen spectrum accumulator
        averageEigenSpectrum += Evals/Evals.sum();
        ++numberOfPointsOnThisPE;
      }  // end if norm > 0
    }  // end GeometryIterator_ loop

//  Compute average spectrum on this PE
    averageEigenSpectrum *= 1/static_cast<double>(numberOfPointsOnThisPE);

//  Compute average eigen spectrum accross all PEs
    std::vector<int> totalNumberOfPoints(1);
    totalNumberOfPoints[0] = numberOfPointsOnThisPE;
    mpiComm.allReduceInPlace(&totalNumberOfPoints[0], 1, eckit::mpi::sum());

    if (numberOfPointsOnThisPE > 0) {
      averageEigenSpectrum = averageEigenSpectrum*
                             static_cast<double>(numberOfPointsOnThisPE)/
                             static_cast<double>(totalNumberOfPoints[0]);
    } else {
      averageEigenSpectrum *= 0;
    }
    mpiComm.allReduceInPlace(averageEigenSpectrum.data(), averageEigenSpectrum.size(),
                             eckit::mpi::sum());
    double averageEigenSpectrumSum = averageEigenSpectrum.sum();
    ASSERT(std::abs(averageEigenSpectrumSum-1.0) < 1e-6);  // by construction trace(C)/nv == 1

//  Convert averageEigenSpectrum to variance explained
    varianceExplained(0) = averageEigenSpectrum(0);
    for (int ii=1; ii < averageEigenSpectrum.size(); ++ii) {
      varianceExplained(ii) = varianceExplained(ii-1) + averageEigenSpectrum(ii);
    }
    oops::Log::test() << "Variance Explained=" << varianceExplained <<std::endl;

//  Compute truncation number based on averageEigenSpectrum
    double frac = 0.0;
    size_t neig = 0;
    for (int ii=0; ii < averageEigenSpectrum.size(); ++ii) {
      frac += averageEigenSpectrum(ii);
      neig += 1;
      if (frac/averageEigenSpectrumSum >= TruncationTollerance) { break; }
    }

//  Compute number of output vectors and retained fraction of the variance
    size_t truncatedNeig = 0;
    if (neig < maxNeigOutput) {
      truncatedNeig = neig;
      frac = frac/averageEigenSpectrumSum;
    } else {
      truncatedNeig = maxNeigOutput;
      frac = averageEigenSpectrum.head(maxNeigOutput).sum()/averageEigenSpectrumSum;
    }
    oops::Log::test() << "Number of retained columns=" << truncatedNeig <<std::endl;
    oops::Log::test() << "Retained fraction of the variance=" << frac <<std::endl;

    return truncatedNeig;
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_SQRTOFVERTLOC_H_
