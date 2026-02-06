/*
 * (C) Crown copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#include <cfloat>
#include <chrono>
#include <utility>

#include "oops/assimilation/ETKFLinearAlgebra.h"

namespace oops {

  Eigen::MatrixXf ETKF_YbRinvYbpI(const Eigen::MatrixXf & Yb,
                                  const Eigen::MatrixXf & YbRinv,
                                  const float scale) {
    // Number of ensembe members in Yb.
    const int nana = Yb.rows();

    // Identity
    const Eigen::VectorXf I = Eigen::VectorXf::Constant(nana, 1.0);

    // Yb R^-1 Yb^T + (nens - 1) I / infl
    const Eigen::MatrixXf YbRinvYbpI =
      YbRinv * Yb.transpose() + I.asDiagonal().toDenseMatrix() * scale;

    return YbRinvYbpI;
  }

  std::tuple<Eigen::VectorXf, Eigen::MatrixXf> Eigendecomposition(const Eigen::MatrixXf & A) {
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(A);
    const Eigen::VectorXf eival = es.eigenvalues().real();
    const Eigen::MatrixXf eivec = es.eigenvectors().real();
    return {eival, eivec};
  }

  std::tuple<Eigen::VectorXf, Eigen::MatrixXf> SingularValueDecomposition(
      const Eigen::MatrixXf & A) {
    const auto svd = Eigen::BDCSVD<Eigen::MatrixXf>(A, Eigen::ComputeFullU);
    const Eigen::VectorXf svals = svd.singularValues().real();
    const Eigen::MatrixXf U = svd.matrixU().real();
    return {svals, U};
  }

  Eigen::MatrixXf ETKF_Pa(const Eigen::VectorXf & eival,
                          const Eigen::MatrixXf & eivec) {
    return eivec * eival.cwiseInverse().asDiagonal() * eivec.transpose();
  }

  Eigen::MatrixXf LETKF_pertWeights(const Eigen::VectorXf & eival,
                                    const Eigen::MatrixXf & eivec) {
    const int nens = eivec.rows();
    return eivec
      * ((nens - 1) * eival.array().inverse()).sqrt().matrix().asDiagonal()
      * eivec.transpose();
  }

  void stoETKF_computeWeights(const Eigen::VectorXf & dy,
                              const Eigen::MatrixXf & YbRinvYbpI,
                              const Eigen::MatrixXf & YbRinv,
                              const Eigen::MatrixXf & YbOrig,
                              const Eigen::SparseMatrix<float> & excludedProjection,
                              const Eigen::SparseMatrix<float> & includedProjection,
                              Eigen::MatrixXf & Wa) {
    // Projecting the excluded members out of matrix YbRinvYbpI
    const Eigen::MatrixXf YbRinvYbpIproj = excludedProjection.transpose()
                                           * YbRinvYbpI * excludedProjection;

    // Eigenvalues and eigenvectors of the above matrix.
    const auto[eival, eivec] = oops::Eigendecomposition(YbRinvYbpIproj);

    // Pa  = [ Yb^T R^-1 Yb + (nens-1)/infl I ] ^-1
    const Eigen::MatrixXf Pa = oops::ETKF_Pa(eival, eivec);

    // Wa = Pa Yb^T R^-1 (dyPert)
    // dyPert = (y_mean + y_pert) - (yb_mean + yb_pert)
    Wa += oops::ETKF_stateWeights<Eigen::MatrixXf>(excludedProjection*Pa
                                                   *excludedProjection.transpose(),
                                                   YbRinv,
                                                   (YbOrig.transpose().colwise()
                                                   + dy)*includedProjection);
  }

  std::tuple<Eigen::MatrixXd, Eigen::VectorXd> ETKF_ensembleIncrement(const Eigen::MatrixXd & Xb,
                                                                      const Eigen::MatrixXd & Wa) {
    const Eigen::MatrixXd Xinc = Xb * Wa;
    return {Xinc, Xinc.rowwise().mean()};
  }

  void ETKF_updateAnalysis(Eigen::MatrixXd & Xa,
                           const Eigen::VectorXd & xinc) {
    Xa.colwise() += xinc;
  }

  Eigen::MatrixXd ETKF_analysisPerturbations(const Eigen::MatrixXd & Xb,
                                             const Eigen::MatrixXd & Xinc,
                                             const Eigen::VectorXd & xinc) {
    return (Xinc + Xb).colwise() - xinc;
  }

  Eigen::VectorXd ETKF_updateEnsembleMean(const Eigen::MatrixXd & Xb,
                                          const Eigen::VectorXd & wa) {
    return Xb * wa;
  }

  Eigen::MatrixXd LETKF_updateEnsemblePerturbation(const Eigen::MatrixXd & Xb,
                                                   const Eigen::MatrixXd & Wa) {
    return Xb * Wa;
  }

  Eigen::MatrixXd GETKF_updateEnsemblePerturbation(const Eigen::MatrixXd & Xb,
                                                   const Eigen::MatrixXd & Xc) {
    return Xb + Xc;
  }

  void ETKF_posteriorInflation(const Eigen::MatrixXd & Xb,
                               Eigen::MatrixXd & Xa,
                               const eckit::Configuration & inflopt) {
    const size_t nens = Xa.cols();

    // RTPP inflation
    const double rtpp = inflopt.getDouble("rtpp", 0.0);
    const double rtps = inflopt.getDouble("rtps", 0.0);
    if (rtpp > 0.0 && rtpp <= 1.0) {
      Xa = (1 - rtpp) * Xa + rtpp * Xb;
    }

    // RTPS inflation
    const double eps = DBL_EPSILON;
    const double rtpsInflMin = 1.0;
    const double rtpsInflMax = 1e30;
    if (rtps > 0.0 && rtps <= 1.0) {
      // posterior spread
      Eigen::ArrayXd asprd = Xa.array().square().rowwise().sum()/(nens - 1);
      asprd = asprd.sqrt();
      asprd = (asprd < eps).select(eps, asprd);  // avoid nan overflow for vars with no spread

      // prior spread
      Eigen::ArrayXd fsprd = Xb.array().square().rowwise().sum()/(nens - 1);
      fsprd = fsprd.sqrt();
      fsprd = (fsprd < eps).select(eps, fsprd);

      // RTPS inflation factor
      Eigen::ArrayXd rtpsInfl = rtps * ((fsprd - asprd)/asprd) + 1;
      rtpsInfl = (rtpsInfl < rtpsInflMin).select(rtpsInflMin, rtpsInfl);
      rtpsInfl = (rtpsInfl > rtpsInflMax).select(rtpsInflMax, rtpsInfl);

      // inflate perturbation matrix
      Xa.array().colwise() *= rtpsInfl;
    }
  }
}  // namespace oops
