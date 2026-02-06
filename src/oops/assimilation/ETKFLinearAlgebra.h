/*
 * (C) Crown copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#ifndef OOPS_ASSIMILATION_ETKFLINEARALGEBRA_H_
#define OOPS_ASSIMILATION_ETKFLINEARALGEBRA_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <tuple>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/generic/VerticalLocEV.h"

namespace oops {

  /// \brief Compute Yb * R^{-1}.
  /// \details
  /// Input:
  ///   Yb: Ensemble perturbations in observation space.
  ///   invVarR: Inverse observation error covariance matrix.
  /// Output:
  ///   Yb * invVarR.
  template <typename OBSERROR>
  Eigen::MatrixXf ETKF_YbRinv(const Eigen::MatrixXf & Yb,
                              const OBSERROR & R) {
    Eigen::MatrixXd YbRinv =  R.localInverseMultiply(Yb);
    return YbRinv.cast<float>();
  }

  /// \brief Compute Yb * R^{-1} * Yb^T + (nens - 1) * I / infl.
  /// \details
  /// Input:
  ///   Yb: Ensemble perturbations in observation space.
  ///   YbRinv: Precomputed Yb * R^{-1}.
  ///   scale: (nens - 1) / infl, where nens is the number of ensemble members
  ///          and infl is a user-defined inflation parameter.
  /// Output:
  ///   Yb * Rinv * Yb^T + scale * I, where I is the identity matrix.
  Eigen::MatrixXf ETKF_YbRinvYbpI(const Eigen::MatrixXf &,
                                  const Eigen::MatrixXf &,
                                  const float);

  /// \brief Perform eigendecomposition on a self-adjoint matrix.
  /// \details
  /// Input:
  ///   A: self-adjoint matrix to be decomposed.
  /// Output:
  ///   {eival, eivec}: tuple of real Eigenvalues and Eigenvectors.
  std::tuple<Eigen::VectorXf, Eigen::MatrixXf> Eigendecomposition(const Eigen::MatrixXf &);

  /// \brief Perform full singular value decomposition using Bidiagonal
  /// Divide and Conquer.
  /// \details
  /// Input:
  ///   A: matrix to be decomposed.
  /// Output:
  ///   {svals, U}: tuple of singular values and left-singular vectors.
  std::tuple<Eigen::VectorXf, Eigen::MatrixXf> SingularValueDecomposition(const Eigen::MatrixXf &);

  /// \brief Compute Pa.
  /// \details
  /// Input:
  ///   eival: Eigenvalues of Yb * R^{-1} * Yb^T + (nens - 1) * I / infl.
  ///   eivec: Eigenvectors of Yb * R^{-1} * Yb^T + (nens - 1) * I / infl.
  /// Output:
  ///   Pa = eivec * eival^{-1/2} * eivec.
  Eigen::MatrixXf ETKF_Pa(const Eigen::VectorXf &,
                          const Eigen::MatrixXf &);

  /// \brief Compute ETKF weights applied to a model state field.
  /// \details
  /// Input:
  ///   Pa: Precomputed Pa matrix.
  ///   YbRinv: Precomputed Yb * R^{-1}.
  ///   y: State in observation space.
  /// Output:
  ///   Pa * Yb * R^{-1} * y
  template <typename T>
  T ETKF_stateWeights(const Eigen::MatrixXf & Pa,
                      const Eigen::MatrixXf & YbRinv,
                      const T & y) {
    return Pa * (YbRinv * y);
  }

  /// \brief Compute LETKF weights applied to a model perturbation field.
  /// \details
  /// Input:
  ///   eival: Eigenvalues of Yb * R^{-1} * Yb^T + (nens - 1) * I / infl.
  ///   eivec: Eigenvectors of Yb * R^{-1} * Yb^T + (nens - 1) * I / infl.
  /// Output:
  ///   LETKF weights = eivec * ((nens - 1) * eival)^{-1/2} * eivec.
  Eigen::MatrixXf LETKF_pertWeights(const Eigen::VectorXf &,
                                    const Eigen::MatrixXf &);

  /// \brief Compute GETKF weights applied to a model perturbation field.
  /// \details
  /// Input:
  ///   eival: Eigenvalues of Yb * R^{-1} * Yb^T + (nens - 1) * I / infl.
  ///   eivec: Eigenvectors of Yb * R^{-1} * Yb^T + (nens - 1) * I / infl.
  ///   Yb: Modulated ensemble perturbations in observation space.
  ///   YbOrig: Original ensemble perturbations in observation space.
  ///   R: Observation error covariance matrix.
  ///   infl: a user-defined inflation parameter.
  /// Output:
  ///   GETKF weights (see code for formula).
  template <typename OBSERROR>
  Eigen::MatrixXf GETKF_pertWeights(const Eigen::VectorXf & eival,
                                    const Eigen::MatrixXf & eivec,
                                    const Eigen::MatrixXf & Yb,
                                    const Eigen::MatrixXf & YbOrig,
                                    const OBSERROR & R,
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
    const Eigen::MatrixXf YbRinv = oops::ETKF_YbRinv(Yb, R);

    // Yb R^-1 YbOrig^T
    const Eigen::MatrixXf YbRinvYbOrig = YbRinv * YbOrig.transpose();

    return -scaledCCT * YbRinvYbOrig;
  }

  /// \brief Compute GETKF weights applied to a model perturbation field
  ///        for a given set of projection matrices.
  /// \details
  /// Input:
  ///   YbRinvYbpI: Matrix equal to Yb * R^{-1} * Yb^T + (nens - 1) * I / infl.
  ///   Yb: Modulated ensemble perturbations in observation space.
  ///   YbOrig: Original ensemble perturbations in observation space.
  ///   R: Observation error covariance matrix.
  ///   excludedProjection: Projection matrix for the excluded
  ///                       subensemble members for cross validation.
  ///   includedProjection: Projection matrix for the included
  ///                       subensemble members for cross validation.
  ///   infl: a user-defined inflation parameter.
  /// Output:
  ///   GETKF weights (see code for formula).
  ///   timings: Vector of timings of internal routines.
  template <typename OBSERROR>
  Eigen::MatrixXf GETKF_pertWeights(const Eigen::MatrixXf & YbRinvYbpI,
                                    const Eigen::MatrixXf & Yb,
                                    const Eigen::MatrixXf & YbOrig,
                                    const OBSERROR & R,
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
    const Eigen::MatrixXf YbRinv = oops::ETKF_YbRinv(Yb, R);

    // Yb R^-1 YbOrig^T
    const Eigen::MatrixXf YbRinvYbOrig = YbRinv * YbOrig.transpose();

    timings = {tE0, tE1, tE2};

    return -excludedProjection * scaledCCT * excludedProjection.transpose()
            * YbRinvYbOrig * includedProjection;
  }

  /// \brief Compute state and perturbation weights for determinstic LETKF.
  /// \details
  /// Input:
  ///   dy: Observation departures.
  ///   Yb: Ensemble perturbations in observation space.
  ///   R: Observation error covariance matrix.
  ///   scale: (nens - 1) / infl, where nens is the number of ensemble members
  ///          and infl is a user-defined inflation parameter.
  /// Output:
  ///    wa [passed by reference]: state weights calculated using ETKF_stateWeights.
  ///    Wa [passed by reference]: perturbation weights calculated using LETKF_pertWeights.
  ///    A vector of times that can be used to profile the performance of this routine.
  template <typename OBSERROR>
  std::vector<std::chrono::time_point<std::chrono::system_clock>>
  detLETKF_computeWeights(const Eigen::VectorXf & dy,
                          const Eigen::MatrixXf & Yb,
                          const OBSERROR & R,
                          const float scale,
                          const bool svd,
                          Eigen::VectorXf & wa,
                          Eigen::MatrixXf & Wa) {
    oops::Log::info() << "!! oops::detLETKF_computeWeights" << std::endl;
    const auto tE0 = std::chrono::system_clock::now();

    // Yb R^-1
    oops::Log::info() << "Calling ETKF_YbRinv" << std::endl;
    const Eigen::MatrixXf YbRinv = oops::ETKF_YbRinv(Yb, R);

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

  /// \brief Compute state and perturbation weights for determinstic GETKF.
  /// \details
  /// Input:
  ///   dy: Observation departures.
  ///   Yb: Modulated ensemble perturbations in observation space.
  ///   YbOrig: Original ensemble perturbations in observation space.
  ///   R: Observation error covariance matrix.
  ///   infl: a user-defined inflation parameter.
  /// Output:
  ///    wa [passed by reference]: state weights calculated using ETKF_stateWeights.
  ///    Wa [passed by reference]: perturbation weights calculated using GETKF_pertWeights.
  ///    A vector of times that can be used to profile the performance of this routine.
  template <typename OBSERROR>
  std::vector<std::chrono::time_point<std::chrono::system_clock>>
  detGETKF_computeWeights(const Eigen::VectorXf & dy,
                          const Eigen::MatrixXf & Yb,
                          const Eigen::MatrixXf & YbOrig,
                          const OBSERROR & R,
                          const float infl,
                          const bool svd,
                          Eigen::VectorXf & wa,
                          Eigen::MatrixXf & Wa) {
    const int nens = YbOrig.rows();

    const auto tE0 = std::chrono::system_clock::now();

    // Yb R^-1
    const Eigen::MatrixXf YbRinv = oops::ETKF_YbRinv(Yb, R);

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
                                 R, infl);

    const auto tE6 = std::chrono::system_clock::now();

    return {tE0, tE1, tE2, tE3, tE4, tE5, tE6};
  }
  /// \brief Computes state and perturbation weights for determinstic ETKF
  ///        for a given set of projection matrices.
  /// \details
  /// Input:
  ///   dy:                        Observation departures.
  ///   Yb:                        Modulated ensemble perturbations in observation space.
  ///   YbRinvYbpI:                Matrix equal to Y^T R^-1 Y + (nens-1)/infl I
  ///   YbRinv:                    Matrix equal to Y^T R^-1
  ///   YbOrig:                    Original ensemble perturbations in observation space.
  ///   R:                         Observation error covariance matrix.
  ///   infl:                      User-defined inflation parameter.
  ///   excludedProjection:        Projection matrix for the excluded
  ///                              subensemble members for cross validation.
  ///   includedProjection:        Projection matrix for the included
  ///                              subensemble members for cross validation.
  ///   computeMeanWeights:        Boolean to compute the mean weights matrix wa.
  /// Output:
  ///    wa [passed by reference]: State weights calculated using ETKF_stateWeights.
  ///    Wa [passed by reference]: Perturbation weights calculated using GETKF_pertWeights.
  ///    A vector of times that can be used to profile the performance of this routine.
  template <typename OBSERROR>
  std::vector<std::chrono::time_point<std::chrono::system_clock>>
  detGETKF_computeWeights(const Eigen::VectorXf & dy,
                          const Eigen::MatrixXf & Yb,
                          const Eigen::MatrixXf & YbRinvYbpI,
                          const Eigen::MatrixXf & YbRinv,
                          const Eigen::MatrixXf & YbOrig,
                          const OBSERROR & R,
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
                                    R,
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

  /// \brief Compute state and perturbation weights for stochastic ETKF.
  /// \details
  /// Input:
  ///   dy: Observation departures.
  ///   Yb: Modulated ensemble perturbations in observation space.
  ///   YbOrig: Original ensemble perturbations in observation space.
  ///   R: Observation error covariance matrix.
  ///   infl: a user-defined inflation parameter.
  ///   svd: if true, perform eigendecomposition of inverse analysis error
  ///   covariance with singular value decomposition algorithm, otherwise use
  ///   eigendecomposition.
  /// Output:
  ///    Wa [passed by reference]: perturbation weights calculated using ETKF_stateWeights.
  template <typename OBSERROR>
  void stoETKF_computeWeights(const Eigen::VectorXf & dy,
                              const Eigen::MatrixXf & Yb,
                              const Eigen::MatrixXf & YbOrig,
                              const OBSERROR & R,
                              const float infl,
                              const bool svd,
                              Eigen::MatrixXf & Wa) {
      const int nens = YbOrig.rows();
      const float scale = (nens - 1) / infl;

      // Yb R^-1
      const Eigen::MatrixXf YbRinv = oops::ETKF_YbRinv(Yb, R);

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
                                                    YbOrig.transpose().colwise()+dy);
  }

  /// \brief Computes state and perturbation weights for stochastic ETKF
  ///        for a given set of projection matrices.
  /// \details
  /// Input:
  ///   dy:                        Observation departures.
  ///   YbRinvYbpI:                Matrix equal to Y^T R^-1 Y + (nens-1)/infl I
  ///   YbRinv:                    Matrix equal to Y^T R^-1
  ///   YbOrig:                    Original ensemble perturbations in observation space.
  ///   excludedProjection:        Projection matrix for the excluded
  ///                              subensemble members for cross validation.
  ///   includedProjection:        Projection matrix for the included
  ///                              subensemble members for cross validation.
  /// Output:
  ///    Wa [passed by reference]: Perturbation weights calculated using ETKF_stateWeights.
  void stoETKF_computeWeights(const Eigen::VectorXf &,
                              const Eigen::MatrixXf &,
                              const Eigen::MatrixXf &,
                              const Eigen::MatrixXf &,
                              const Eigen::SparseMatrix<float> &,
                              const Eigen::SparseMatrix<float> &,
                              Eigen::MatrixXf &);

  /// \brief Compute ensemble increment.
  /// \details
  /// Input:
  ///   Xb: Ensemble perturbations in model space.
  ///   Wa: Precomputed perturbation weights.
  /// Output:
  ///   Xinc: Xb * Wa.
  ///   xinc: row-wise mean of Xinc.
  std::tuple<Eigen::MatrixXd, Eigen::VectorXd> ETKF_ensembleIncrement(const Eigen::MatrixXd &,
                                                                      const Eigen::MatrixXd &);

  /// \brief Update analysis state.
  /// \details
  /// Input:
  ///   xinc: Precomputed row-wise mean of (Xb * Wa)
  /// Output:
  ///   Xa [passed by reference]: State Xa multiplied column-wise by xinc.
  void ETKF_updateAnalysis(Eigen::MatrixXd &, const Eigen::VectorXd &);

  /// \brief Compute analysis perturbations.
  /// \details
  /// Input:
  ///   Xb: Ensemble perturbations in model space.
  ///   Xinc: Precomputed Xb * Wa.
  ///   xinc: Precomputed row-wise mean of Xinc.
  /// Output:
  ///   Xa: Analysis in model space (see code for formula).
  Eigen::MatrixXd ETKF_analysisPerturbations(const Eigen::MatrixXd &,
                                             const Eigen::MatrixXd &,
                                             const Eigen::VectorXd &);

  /// \brief Update ensemble mean.
  /// \details
  /// Input:
  ///   Xb: Ensemble perturbations in model space.
  ///   wa: Precomputed state weights.
  /// Output:
  ///   Xb * wa.
  Eigen::VectorXd ETKF_updateEnsembleMean(const Eigen::MatrixXd &,
                                          const Eigen::VectorXd &);

  /// \brief Update LETKF ensemble perturbations.
  /// \details
  /// Input:
  ///   Xb: Ensemble perturbations in model space.
  ///   Wa: Precomputed perturbation weights.
  /// Output:
  ///   Xb * Wa.
  Eigen::MatrixXd LETKF_updateEnsemblePerturbation(const Eigen::MatrixXd &,
                                                   const Eigen::MatrixXd &);

  /// \brief Update GETKF ensemble perturbations.
  /// \details
  /// Input:
  ///   Xb: Ensemble perturbations in model space.
  ///   Xc: Perturbation matrix.
  /// Output:
  ///   Xb + Xc.
  Eigen::MatrixXd GETKF_updateEnsemblePerturbation(const Eigen::MatrixXd &,
                                                   const Eigen::MatrixXd &);

  /// \brief Apply posterior inflation (RTPP or RTPS) to a state.
  /// \details
  /// Input:
  ///   Xb: Ensemble perturbations in model space.
  ///   inflopt: inflation configuration options (see code for details).
  /// Output:
  ///   Xa [passed by reference]: Output state in model space.
  void ETKF_posteriorInflation(const Eigen::MatrixXd &,
                               Eigen::MatrixXd &,
                               const eckit::Configuration &);

  /// \brief Perform stochastic ensemble transform.
  /// \details
  /// Input:
  ///   bkg_pert: Ensemble background perturbations.
  ///   geomIter: Model geometry iterator.
  ///   Wa: Precomputed perturbation weights.
  ///   inflopt: inflation configuration options (see code for details).
  ///   vertLoc [optional]: vertical localization; triggers ensemble modulation
  ///                       if passed in.
  /// Output:
  ///   ana_pert [passed by reference]: Ensemble analysis perturbations.
  template <typename MODEL, typename IncSet, typename GeomIt>
  void stoETKF_applyWeights(const IncSet & bkg_pert,
                            IncSet & ana_pert,
                            const GeomIt & geomIter,
                            const Eigen::MatrixXf & Wa,
                            const eckit::Configuration & inflopt,
                            const VerticalLocEV<MODEL> * const vertLoc = nullptr) {
    // Loop through analysis times.
    for (size_t itime = 0; itime < bkg_pert.time_size(); ++itime) {
      // Original Xb.
      Eigen::MatrixXd Xb;
      bkg_pert.packEigen(Xb, geomIter, itime);

      // Ensemble perturbations in model space.
      // This is either the original or the modulated Xb.
      const Eigen::MatrixXd Xc = vertLoc ?
        vertLoc->modulateIncrement(bkg_pert, geomIter, itime) : Xb;

      // Compute the ensemnble increment using either original or modulated Xb.
      // Xinc = Xb * Wa and xinc is the row-wise mean of Xinc.
      const auto[Xinc, xinc] =
        oops::ETKF_ensembleIncrement(Xc, Wa.cast<double>());

      // Generate analysis perturbations for inflation.
      Eigen::MatrixXd Xa = oops::ETKF_analysisPerturbations(Xb, Xinc, xinc);

      // Apply posterior inflation.
      oops::ETKF_posteriorInflation(Xb, Xa, inflopt);

      // Update analysis.
      oops::ETKF_updateAnalysis(Xa, xinc);

      // Update analysis perturbations.
      ana_pert.setEigen(Xa, geomIter, itime);
    }
  }

  /// \brief Perform deterministic ensemble transform.
  /// \details
  /// Input:
  ///   bkg_pert: Ensemble background perturbations.
  ///   geomIter: Model geometry iterator.
  ///   wa: Precomputed state weights.
  ///   Wa: Precomputed perturbation weights.
  ///   inflopt: inflation configuration options (see code for details).
  ///   vertLoc [optional]: vertical localization; triggers ensemble modulation
  ///                       if passed in.
  /// Output:
  ///   ana_pert [passed by reference]: Ensemble analysis perturbations.
  template <typename MODEL, typename IncSet, typename GeomIt>
  void detETKF_applyWeights(const IncSet & bkg_pert,
                            IncSet & ana_pert,
                            const GeomIt & geomIter,
                            const Eigen::VectorXf & wa,
                            const Eigen::MatrixXf & Wa,
                            const eckit::Configuration & inflopt,
                            const VerticalLocEV<MODEL> * const vertLoc = nullptr) {
    // Loop through analysis times.
    for (size_t itime = 0; itime < bkg_pert.time_size(); ++itime) {
      // Original Xb.
      Eigen::MatrixXd Xb;
      bkg_pert.packEigen(Xb, geomIter, itime);

      // Ensemble perturbations in model space.
      // This is either the original or the modulated Xb.
      const Eigen::MatrixXd Xc = vertLoc ?
        vertLoc->modulateIncrement(bkg_pert, geomIter, itime) : Xb;

      // Update ensemble mean using either original or modulated Xb.
      const Eigen::VectorXd xa = oops::ETKF_updateEnsembleMean(Xc, wa.cast<double>());

      // Update ensemble perturbations using either original or modulated Xb.
      Eigen::MatrixXd Xa = vertLoc ?
        oops::GETKF_updateEnsemblePerturbation(Xb, Xc * Wa.cast<double>()) :
        oops::LETKF_updateEnsemblePerturbation(Xb, Wa.cast<double>());

      // Apply posterior inflation.
      oops::ETKF_posteriorInflation(Xb, Xa, inflopt);

      // Update analysis.
      oops::ETKF_updateAnalysis(Xa, xa);

      // Update analysis perturbations.
      ana_pert.setEigen(Xa, geomIter, itime);
    }
  }

}  // namespace oops

#endif  // OOPS_ASSIMILATION_ETKFLINEARALGEBRA_H_
