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
  Eigen::MatrixXd ETKF_YbRinv(const Eigen::MatrixXf &,
                              const Eigen::VectorXd &);

  /// \brief Compute Yb * R^{-1} * Yb^T + (nens - 1) * I / infl.
  /// \details
  /// Input:
  ///   Yb: Ensemble perturbations in observation space.
  ///   YbRinv: Precomputed Yb * R^{-1}.
  ///   scale: (nens - 1) / infl, where nens is the number of ensemble members
  ///          and infl is a user-defined inflation parameter.
  /// Output:
  ///   Yb * Rinv * Yb^T + scale * I, where I is the identity matrix.
  Eigen::MatrixXd ETKF_YbRinvYbpI(const Eigen::MatrixXf &,
                                  const Eigen::MatrixXd &,
                                  const double);

  /// \brief Perform Eigendecomposition.
  /// \details
  /// Input:
  ///   A: matrix to be decomposed.
  /// Output:
  ///   {eival, eivec}: tuple of Eigenvalues and Eigenvectors.
  std::tuple<Eigen::VectorXd, Eigen::MatrixXd> Eigendecomposition(const Eigen::MatrixXd &);

  /// \brief Compute Pa.
  /// \details
  /// Input:
  ///   eival: Eigenvalues of Yb * R^{-1} * Yb^T + (nens - 1) * I / infl.
  ///   eivec: Eigenvectors of Yb * R^{-1} * Yb^T + (nens - 1) * I / infl.
  /// Output:
  ///   Pa = eivec * eival^{-1/2} * eivec.
  Eigen::MatrixXd ETKF_Pa(const Eigen::VectorXd &,
                          const Eigen::MatrixXd &);

  /// \brief Compute ETKF weights applied to a model state field.
  /// \details
  /// Input:
  ///   Pa: Precomputed Pa matrix.
  ///   YbRinv: Precomputed Yb * R^{-1}.
  ///   y: State in observation space.
  /// Output:
  ///   Pa * Yb * R^{-1} * y
  template <typename T>
  T ETKF_stateWeights(const Eigen::MatrixXd & Pa,
                      const Eigen::MatrixXd & YbRinv,
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
  Eigen::MatrixXd LETKF_pertWeights(const Eigen::VectorXd &,
                                    const Eigen::MatrixXd &);

  /// \brief Compute GETKF weights applied to a model perturbation field.
  /// \details
  /// Input:
  ///   eival: Eigenvalues of Yb * R^{-1} * Yb^T + (nens - 1) * I / infl.
  ///   eivec: Eigenvectors of Yb * R^{-1} * Yb^T + (nens - 1) * I / infl.
  ///   Yb: Modulated ensemble perturbations in observation space.
  ///   YbOrig: Original ensemble perturbations in observation space.
  ///   invVarR: Inverse observation error covariance matrix.
  ///   infl: a user-defined inflation parameter.
  /// Output:
  ///   GETKF weights (see code for formula).
  Eigen::MatrixXd GETKF_pertWeights(const Eigen::VectorXd &,
                                    const Eigen::MatrixXd &,
                                    const Eigen::MatrixXf &,
                                    const Eigen::MatrixXf &,
                                    const Eigen::VectorXd &,
                                    const float);

  /// \brief Compute state and perturbation weights for determinstic LETKF.
  /// \details
  /// Input:
  ///   dy: Observation departures.
  ///   Yb: Ensemble perturbations in observation space.
  ///   invVarR: Inverse observation error covariance matrix.
  ///   scale: (nens - 1) / infl, where nens is the number of ensemble members
  ///          and infl is a user-defined inflation parameter.
  /// Output:
  ///    wa [passed by reference]: state weights calculated using ETKF_stateWeights.
  ///    Wa [passed by reference]: perturbation weights calculated using LETKF_pertWeights.
  ///    A vector of times that can be used to profile the performance of this routine.
  std::vector<std::chrono::time_point<std::chrono::system_clock>>
  detLETKF_computeWeights(const Eigen::VectorXd &,
                          const Eigen::MatrixXf &,
                          const Eigen::VectorXd &,
                          const double,
                          Eigen::VectorXd &,
                          Eigen::MatrixXd &);

  /// \brief Compute state and perturbation weights for determinstic GETKF.
  /// \details
  /// Input:
  ///   dy: Observation departures.
  ///   Yb: Modulated ensemble perturbations in observation space.
  ///   YbOrig: Original ensemble perturbations in observation space.
  ///   invVarR: Inverse observation error covariance matrix.
  ///   infl: a user-defined inflation parameter.
  /// Output:
  ///    wa [passed by reference]: state weights calculated using ETKF_stateWeights.
  ///    Wa [passed by reference]: perturbation weights calculated using GETKF_pertWeights.
  ///    A vector of times that can be used to profile the performance of this routine.
  std::vector<std::chrono::time_point<std::chrono::system_clock>>
  detGETKF_computeWeights(const Eigen::VectorXd &,
                          const Eigen::MatrixXf &,
                          const Eigen::MatrixXf &,
                          const Eigen::VectorXd &,
                          const double,
                          Eigen::VectorXd &,
                          Eigen::MatrixXd &);

  /// \brief Compute state and perturbation weights for stochastic ETKF.
  /// \details
  /// Input:
  ///   dy: Observation departures.
  ///   Yb: Modulated ensemble perturbations in observation space.
  ///   YbOrig: Original ensemble perturbations in observation space.
  ///   invVarR: Inverse observation error covariance matrix.
  ///   infl: a user-defined inflation parameter.
  /// Output:
  ///    Wa [passed by reference]: perturbation weights calculated using ETKF_stateWeights.
  void stoETKF_computeWeights(const Eigen::VectorXd &,
                              const Eigen::MatrixXf &,
                              const Eigen::MatrixXf &,
                              const Eigen::VectorXd &,
                              const double,
                              Eigen::MatrixXd &);

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
                            const Eigen::MatrixXd & Wa,
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
        oops::ETKF_ensembleIncrement(Xc, Wa);

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
                            const Eigen::VectorXd & wa,
                            const Eigen::MatrixXd & Wa,
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
      const Eigen::VectorXd xa = oops::ETKF_updateEnsembleMean(Xc, wa);

      // Update ensemble perturbations using either original or modulated Xb.
      Eigen::MatrixXd Xa = vertLoc ?
        oops::GETKF_updateEnsemblePerturbation(Xb, Xc * Wa) :
        oops::LETKF_updateEnsemblePerturbation(Xb, Wa);

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
