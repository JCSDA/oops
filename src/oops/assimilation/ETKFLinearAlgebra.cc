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

  Eigen::MatrixXf ETKF_YbRinv(const Eigen::MatrixXf & Yb,
                              const Eigen::VectorXf & invVarR) {
    return Yb * invVarR.asDiagonal();
  }

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

  Eigen::MatrixXf GETKF_pertWeights(const Eigen::VectorXf & eival,
                                    const Eigen::MatrixXf & eivec,
                                    const Eigen::MatrixXf & Yb,
                                    const Eigen::MatrixXf & YbOrig,
                                    const Eigen::VectorXf & invVarR,
                                    const float infl) {
    // Normalisation
    const int nens = YbOrig.rows();
    const float norm = 1.0 / (nens - 1.0);

    // Identity
    const int nana = Yb.rows();
    const Eigen::VectorXf I = Eigen::VectorXf::Constant(nana, 1.0);

    // Account for division by zero than occurs when an eigenvalue
    // is exactly equal to 1 / (infl * norm).
    const Eigen::VectorXf eivalmI = (eival - I / (infl * norm)).cwiseAbs();
    const Eigen::VectorXf epsilon = Eigen::VectorXf::Constant(nana, 1.0e-6);
    const Eigen::VectorXf eivalmIsafe =
      (eivalmI.array() < epsilon.array()).select(epsilon, eivalmI);

    // (I - Gamma^{-1/2} / (nens - 1)) * (Gamma - (nens - 1) I / rho)
    const Eigen::VectorXf diag = (I - (norm * eival).cwiseInverse().cwiseAbs().cwiseSqrt()).
      cwiseProduct(eivalmIsafe.cwiseInverse());

    // C ((I - Gamma^{-1/2} / (nens - 1)) * (Gamma - (nens - 1) I / rho)) C^T
    const Eigen::MatrixXf scaledCCT = eivec * diag.asDiagonal() * eivec.transpose();

    // Yb R^-1
    const Eigen::MatrixXf YbRinv = oops::ETKF_YbRinv(Yb, invVarR);

    // Yb R^-1 YbOrig^T
    const Eigen::MatrixXf YbRinvYbOrig = YbRinv * YbOrig.transpose();

    return -scaledCCT * YbRinvYbOrig;
  }

  Eigen::MatrixXf GETKF_pertWeights(const Eigen::MatrixXf & YbRinvYbpI,
                                    const Eigen::MatrixXf & Yb,
                                    const Eigen::MatrixXf & YbOrig,
                                    const Eigen::VectorXf & invVarR,
                                    const Eigen::SparseMatrix<float> & excludedProjection,
                                    const Eigen::SparseMatrix<float> & includedProjection,
                                    const float infl,
                                    std::vector<std::chrono::time_point
                                    <std::chrono::system_clock>> & timings) {
    const auto tE0 = std::chrono::system_clock::now();

    // Projecting the excluded members out of the A matrix
    const Eigen::MatrixXf YbRinvYbpIproj = excludedProjection.transpose()
                                           * YbRinvYbpI * excludedProjection;

    const auto tE1 = std::chrono::system_clock::now();

    // Eigenvalues and eigenvectors of the projected A matrix
    const auto[eival, eivec] = oops::Eigendecomposition(YbRinvYbpIproj);

    const auto tE2 = std::chrono::system_clock::now();

    // Normalisation
    const int nhat = excludedProjection.cols();
    const float norm = 1.0 / (nhat - 1.0);

    // Identity
    const Eigen::VectorXf I = Eigen::VectorXf::Constant(nhat, 1.0);

    // Account for division by zero than occurs when an eigenvalue
    // is exactly equal to 1 / (infl * norm).
    const Eigen::VectorXf eivalmI = (eival - I / (infl * norm)).cwiseAbs();
    const Eigen::VectorXf epsilon = Eigen::VectorXf::Constant(nhat, 1.0e-6);
    const Eigen::VectorXf eivalmIsafe =
      (eivalmI.array() < epsilon.array()).select(epsilon, eivalmI);

    // (I - Gamma^{-1/2} / (nens - 1)) * (Gamma - (nens - 1) I / rho)
    const Eigen::VectorXf diag = (I - (norm * eival).cwiseInverse().cwiseAbs().cwiseSqrt()).
      cwiseProduct(eivalmIsafe.cwiseInverse());

    // C ((I - Gamma^{-1/2} / (nens - 1)) * (Gamma - (nens - 1) I / rho)) C^T
    const Eigen::MatrixXf scaledCCT = eivec * diag.asDiagonal() * eivec.transpose();

    // Yb R^-1
    const Eigen::MatrixXf YbRinv = oops::ETKF_YbRinv(Yb, invVarR);

    // Yb R^-1 YbOrig^T
    const Eigen::MatrixXf YbRinvYbOrig = YbRinv * YbOrig.transpose();

    timings = {tE0, tE1, tE2};

    return -excludedProjection * scaledCCT * excludedProjection.transpose()
            * YbRinvYbOrig * includedProjection;
  }

  std::vector<std::chrono::time_point<std::chrono::system_clock>>
  detLETKF_computeWeights(const Eigen::VectorXf & dy,
                          const Eigen::MatrixXf & Yb,
                          const Eigen::VectorXf & invVarR,
                          const float scale,
                          const bool svd,
                          Eigen::VectorXf & wa,
                          Eigen::MatrixXf & Wa) {
    const auto tE0 = std::chrono::system_clock::now();

    // Yb R^-1
    const Eigen::MatrixXf YbRinv = oops::ETKF_YbRinv(Yb, invVarR);

    // YbRinvYbp = Y^T R^-1 Y + (nens-1)/infl I
    const Eigen::MatrixXf YbRinvYbpI = oops::ETKF_YbRinvYbpI(Yb, YbRinv, scale);

    const auto tE1 = std::chrono::system_clock::now();

    // Eigenvalues and eigenvectors of the above matrix, optionally computed with
    // SVD (matrix should be real symmetric positive definite so these are equivalent).
    const auto[eival, eivec] =
        svd ? oops::SingularValueDecomposition(YbRinvYbpI) : oops::Eigendecomposition(YbRinvYbpI);

    const auto tE2 = std::chrono::system_clock::now();

    // Pa = [Yb R^-1 Yb^T + (nens - 1)/infl I]^-1
    const Eigen::MatrixXf Pa = oops::ETKF_Pa(eival, eivec);

    const auto tE3 = std::chrono::system_clock::now();

    // wa = Pa Yb R^-1 dy
    wa = oops::ETKF_stateWeights<Eigen::VectorXf>(Pa, YbRinv, dy);

    // Wa = sqrt(nens - 1) * eivec * eival^{1/2} * eivec^T
    Wa = oops::LETKF_pertWeights(eival, eivec);

    const auto tE4 = std::chrono::system_clock::now();

    return {tE0, tE1, tE2, tE3, tE4};
  }

  std::vector<std::chrono::time_point<std::chrono::system_clock>>
  detGETKF_computeWeights(const Eigen::VectorXf & dy,
                          const Eigen::MatrixXf & Yb,
                          const Eigen::MatrixXf & YbOrig,
                          const Eigen::VectorXf & invVarR,
                          const float infl,
                          const bool svd,
                          Eigen::VectorXf & wa,
                          Eigen::MatrixXf & Wa) {
    const int nens = YbOrig.rows();

    const auto tE0 = std::chrono::system_clock::now();

    // Yb R^-1
    const Eigen::MatrixXf YbRinv = oops::ETKF_YbRinv(Yb, invVarR);

    const auto tE1 = std::chrono::system_clock::now();

    const float scale = (nens - 1) / infl;

    // Yb R^-1 Yb^T + (nens - 1) I / infl
    const Eigen::MatrixXf YbRinvYbpI = oops::ETKF_YbRinvYbpI(Yb, YbRinv, scale);

    const auto tE2 = std::chrono::system_clock::now();

    // Eigenvalues and eigenvectors of the above matrix, optionally computed with
    // SVD (matrix should be real symmetric positive definite so these are equivalent).
    const auto[eival, eivec] =
        svd ? oops::SingularValueDecomposition(YbRinvYbpI) : oops::Eigendecomposition(YbRinvYbpI);

    const auto tE3 = std::chrono::system_clock::now();

    // Pa = [Yb R^-1 Yb^T + (nens - 1)/infl I]^-1
    const Eigen::MatrixXf Pa = oops::ETKF_Pa(eival, eivec);

    const auto tE4 = std::chrono::system_clock::now();

    // wa = Pa Yb R^-1 dy
    wa = oops::ETKF_stateWeights<Eigen::VectorXf>(Pa, YbRinv, dy);

    const auto tE5 = std::chrono::system_clock::now();

    // Wa = eivec ((I - eival^{-1/2} / (nens - 1)) * (eival - (nens - 1) I / rho)) eivec^T
    //      * Yb R^-1 YbOrig^T
    Wa = oops::GETKF_pertWeights(eival, eivec,
                                 Yb, YbOrig,
                                 invVarR, infl);

    const auto tE6 = std::chrono::system_clock::now();

    return {tE0, tE1, tE2, tE3, tE4, tE5, tE6};
  }

  std::vector<std::chrono::time_point<std::chrono::system_clock>>
  detGETKF_computeWeights(const Eigen::VectorXf & dy,
                          const Eigen::MatrixXf & Yb,
                          const Eigen::MatrixXf & YbRinvYbpI,
                          const Eigen::MatrixXf & YbRinv,
                          const Eigen::MatrixXf & YbOrig,
                          const Eigen::VectorXf & invVarR,
                          const float infl,
                          const Eigen::SparseMatrix<float> & excludedProjection,
                          const Eigen::SparseMatrix<float> & includedProjection,
                          const bool computeMeanWeights,
                          Eigen::VectorXf & wa,
                          Eigen::MatrixXf & Wa) {
    const auto tE0 = std::chrono::system_clock::now();

    if (computeMeanWeights) {
      // Eigenvalues and eigenvectors of the A matrix
      const auto[eival, eivec] = oops::Eigendecomposition(YbRinvYbpI);

      // Pa = [Yb R^-1 Yb^T + (nhat - 1)/infl I]^-1
      const Eigen::MatrixXf Pa = oops::ETKF_Pa(eival, eivec);

      // wa = Pa Yb R^-1 dy
      wa = oops::ETKF_stateWeights<Eigen::VectorXf>(Pa, YbRinv, dy);
    }

    // Wa = eivec ((I - eival^{-1/2} / (nens - 1)) * (eival - (nens - 1) I / rho)) eivec^T
    //      * Yb R^-1 YbOrig^T
    std::vector<std::chrono::time_point<std::chrono::system_clock>> timings;
    Wa += oops::GETKF_pertWeights(YbRinvYbpI,
                                  Yb, YbOrig,
                                  invVarR,
                                  excludedProjection, includedProjection,
                                  infl,
                                  timings);

    const auto tE4 = std::chrono::system_clock::now();

    // Unpacking timings
    const auto tE1 = timings[0];
    const auto tE2 = timings[1];
    const auto tE3 = timings[2];

    return {tE0, tE1, tE2, tE3, tE4};
  }

  void stoETKF_computeWeights(const Eigen::VectorXf & dy,
                              const Eigen::MatrixXf & Yb,
                              const Eigen::MatrixXf & YbOrig,
                              const Eigen::VectorXf & invVarR,
                              const float infl,
                              const bool svd,
                              Eigen::MatrixXf & Wa) {
    const int nens = YbOrig.rows();
    const float scale = (nens - 1) / infl;

    // Yb R^-1
    const Eigen::MatrixXf YbRinv = oops::ETKF_YbRinv(Yb, invVarR);

    // YbRinvYbp = Y^T R^-1 Y + (nens-1)/infl I
    const Eigen::MatrixXf YbRinvYbpI = oops::ETKF_YbRinvYbpI(Yb, YbRinv, scale);

    // Eigenvalues and eigenvectors of the above matrix, optionally computed with
    // SVD (matrix should be real symmetric positive definite so these are equivalent).
    const auto[eival, eivec] =
        svd ? oops::SingularValueDecomposition(YbRinvYbpI) : oops::Eigendecomposition(YbRinvYbpI);

    // Pa  = [ Yb^T R^-1 Yb + (nens-1)/infl I ] ^-1
    const Eigen::MatrixXf Pa = oops::ETKF_Pa(eival, eivec);

    // Wa = Pa Yb^T R^-1 (dyPert)
    // dyPert = (y_mean + y_pert) - (yb_mean + yb_pert)
    Wa = oops::ETKF_stateWeights<Eigen::MatrixXf>(Pa,
                                                  YbRinv,
                                                  YbOrig.transpose().colwise() + dy);
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
