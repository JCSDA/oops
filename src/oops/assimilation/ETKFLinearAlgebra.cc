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

  Eigen::MatrixXd ETKF_YbRinv(const Eigen::MatrixXf & Yb,
                              const Eigen::VectorXd & invVarR) {
    return Yb.cast<double>() * invVarR.asDiagonal();
  }

  Eigen::MatrixXd ETKF_YbRinvYbpI(const Eigen::MatrixXf & Yb,
                                  const Eigen::MatrixXd & YbRinv,
                                  const double scale) {
    // Number of ensembe members in Yb.
    const int nana = Yb.rows();

    // Identity
    const Eigen::VectorXd I = Eigen::VectorXd::Constant(nana, 1.0);

    // Yb R^-1 Yb^T + (nens - 1) I / infl
    const Eigen::MatrixXd YbRinvYbpI =
      YbRinv * Yb.cast<double>().transpose() + I.asDiagonal().toDenseMatrix() * scale;

    return YbRinvYbpI;
  }

  std::tuple<Eigen::VectorXd, Eigen::MatrixXd> Eigendecomposition(const Eigen::MatrixXd & A) {
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
    const Eigen::VectorXd eival = es.eigenvalues().real();
    const Eigen::MatrixXd eivec = es.eigenvectors().real();
    return {eival, eivec};
  }

  Eigen::MatrixXd ETKF_Pa(const Eigen::VectorXd & eival,
                          const Eigen::MatrixXd & eivec) {
    return eivec * eival.cwiseInverse().asDiagonal() * eivec.transpose();
  }

  Eigen::MatrixXd LETKF_pertWeights(const Eigen::VectorXd & eival,
                                    const Eigen::MatrixXd & eivec) {
    const int nens = eivec.rows();
    return eivec
      * ((nens - 1) * eival.array().inverse()).sqrt().matrix().asDiagonal()
      * eivec.transpose();
  }

  Eigen::MatrixXd GETKF_pertWeights(const Eigen::VectorXd & eival,
                                    const Eigen::MatrixXd & eivec,
                                    const Eigen::MatrixXf & Yb,
                                    const Eigen::MatrixXf & YbOrig,
                                    const Eigen::VectorXd & invVarR,
                                    const float infl) {
    // Normalisation
    const int nens = YbOrig.rows();
    const float norm = 1.0 / (nens - 1.0);

    // Identity
    const int nana = Yb.rows();
    const Eigen::VectorXd I = Eigen::VectorXd::Constant(nana, 1.0);

    // Account for division by zero than occurs when an eigenvalue
    // is exactly equal to 1 / (infl * norm).
    const Eigen::VectorXd eivalmI = (eival - I / (infl * norm)).cwiseAbs();
    const Eigen::VectorXd epsilon = Eigen::VectorXd::Constant(nana, 1.0e-6);
    const Eigen::VectorXd eivalmIsafe =
      (eivalmI.array() < epsilon.array()).select(epsilon, eivalmI);

    // (I - Gamma^{-1/2} / (nens - 1)) * (Gamma - (nens - 1) I / rho)
    const Eigen::VectorXd diag = (I - (norm * eival).cwiseInverse().cwiseAbs().cwiseSqrt()).
      cwiseProduct(eivalmIsafe.cwiseInverse());

    // C ((I - Gamma^{-1/2} / (nens - 1)) * (Gamma - (nens - 1) I / rho)) C^T
    const Eigen::MatrixXd scaledCCT = eivec * diag.asDiagonal() * eivec.transpose();

    // Yb R^-1
    const Eigen::MatrixXd YbRinv = oops::ETKF_YbRinv(Yb, invVarR);

    // Yb R^-1 YbOrig^T
    const Eigen::MatrixXd YbRinvYbOrig = YbRinv * YbOrig.cast<double>().transpose();

    return -scaledCCT * YbRinvYbOrig;
  }

  std::vector<std::chrono::time_point<std::chrono::system_clock>>
  detLETKF_computeWeights(const Eigen::VectorXd & dy,
                          const Eigen::MatrixXf & Yb,
                          const Eigen::VectorXd & invVarR,
                          const double scale,
                          Eigen::VectorXd & wa,
                          Eigen::MatrixXd & Wa) {
    const auto tE0 = std::chrono::system_clock::now();

    // Yb R^-1
    const Eigen::MatrixXd YbRinv = oops::ETKF_YbRinv(Yb, invVarR);

    // YbRinvYbp = Y^T R^-1 Y + (nens-1)/infl I
    const Eigen::MatrixXd YbRinvYbpI = oops::ETKF_YbRinvYbpI(Yb, YbRinv, scale);

    const auto tE1 = std::chrono::system_clock::now();

    // Eigenvalues and eigenvectors of the above matrix.
    const auto[eival, eivec] = oops::Eigendecomposition(YbRinvYbpI);

    const auto tE2 = std::chrono::system_clock::now();

    // Pa = [Yb R^-1 Yb^T + (nens - 1)/infl I]^-1
    const Eigen::MatrixXd Pa = oops::ETKF_Pa(eival, eivec);

    const auto tE3 = std::chrono::system_clock::now();

    // wa = Pa Yb R^-1 dy
    wa = oops::ETKF_stateWeights<Eigen::VectorXd>(Pa, YbRinv, dy);

    // Wa = sqrt(nens - 1) * eivec * eival^{1/2} * eivec^T
    Wa = oops::LETKF_pertWeights(eival, eivec);

    const auto tE4 = std::chrono::system_clock::now();

    return {tE0, tE1, tE2, tE3, tE4};
  }

  std::vector<std::chrono::time_point<std::chrono::system_clock>>
  detGETKF_computeWeights(const Eigen::VectorXd & dy,
                          const Eigen::MatrixXf & Yb,
                          const Eigen::MatrixXf & YbOrig,
                          const Eigen::VectorXd & invVarR,
                          const double infl,
                          Eigen::VectorXd & wa,
                          Eigen::MatrixXd & Wa) {
    const int nens = YbOrig.rows();

    const auto tE0 = std::chrono::system_clock::now();

    // Yb R^-1
    const Eigen::MatrixXd YbRinv = oops::ETKF_YbRinv(Yb, invVarR);

    const auto tE1 = std::chrono::system_clock::now();

    const double scale = (nens - 1) / infl;

    // Yb R^-1 Yb^T + (nens - 1) I / infl
    const Eigen::MatrixXd YbRinvYbpI = oops::ETKF_YbRinvYbpI(Yb, YbRinv, scale);

    const auto tE2 = std::chrono::system_clock::now();

    // Eigenvalues and eigenvectors of the above matrix.
    const auto[eival, eivec] = oops::Eigendecomposition(YbRinvYbpI);

    const auto tE3 = std::chrono::system_clock::now();

    // Pa = [Yb R^-1 Yb^T + (nens - 1)/infl I]^-1
    const Eigen::MatrixXd Pa = oops::ETKF_Pa(eival, eivec);

    const auto tE4 = std::chrono::system_clock::now();

    // wa = Pa Yb R^-1 dy
    wa = oops::ETKF_stateWeights<Eigen::VectorXd>(Pa, YbRinv, dy);

    const auto tE5 = std::chrono::system_clock::now();

    // Wa = eivec ((I - eival^{-1/2} / (nens - 1)) * (eival - (nens - 1) I / rho)) eivec^T
    //      * Yb R^-1 YbOrig^T
    Wa = oops::GETKF_pertWeights(eival, eivec,
                                 Yb, YbOrig,
                                 invVarR, infl);

    const auto tE6 = std::chrono::system_clock::now();

    return {tE0, tE1, tE2, tE3, tE4, tE5, tE6};
  }

  void stoETKF_computeWeights(const Eigen::VectorXd & dy,
                              const Eigen::MatrixXf & Yb,
                              const Eigen::MatrixXf & YbOrig,
                              const Eigen::VectorXd & invVarR,
                              const double infl,
                              Eigen::MatrixXd & Wa) {
    const int nens = YbOrig.rows();
    const double scale = (nens - 1) / infl;

    // Yb R^-1
    const Eigen::MatrixXd YbRinv = oops::ETKF_YbRinv(Yb, invVarR);

    // YbRinvYbp = Y^T R^-1 Y + (nens-1)/infl I
    const Eigen::MatrixXd YbRinvYbpI = oops::ETKF_YbRinvYbpI(Yb, YbRinv, scale);

    // Eigenvalues and eigenvectors of the above matrix.
    const auto[eival, eivec] = oops::Eigendecomposition(YbRinvYbpI);

    // Pa  = [ Yb^T R^-1 Yb + (nens-1)/infl I ] ^-1
    const Eigen::MatrixXd Pa = oops::ETKF_Pa(eival, eivec);

    // Wa = Pa Yb^T R^-1 (dyPert)
    // dyPert = (y_mean + y_pert) - (yb_mean + yb_pert)
    Wa = oops::ETKF_stateWeights<Eigen::MatrixXd>(Pa,
                                                  YbRinv,
                                                  YbOrig.cast<double>().transpose().colwise() + dy);
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
