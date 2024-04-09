/*
 * (C) Copyright 2020 UCAR.
 * (C) Crown Copyright 2023, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_MINIMIZERUTILS_H_
#define OOPS_ASSIMILATION_MINIMIZERUTILS_H_

#include <memory>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/HtRinvHMatrix.h"
#include "oops/base/StructuredGridWriter.h"
#include "oops/util/Logger.h"

namespace oops {

/// Prints to Log::info gradient reduction \p grad and normalized gradient reduction \p norm
/// for iteration \p iteration
void printNormReduction(int iteration, const double & grad, const double & norm);
/// Prints to Log::info cost function values for \p costJ, \p costJb, \p costJoJc for
/// iteration \p iteration
void printQuadraticCostFunction(int iteration, const double & costJ, const double & costJb,
                                const double & costJoJc);

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void writeIncrement(const eckit::Configuration & config,
                    const ControlIncrement<MODEL, OBS> & dx, const int & loop) {
// Write out the increment

  if (config.has("online diagnostics")) {
    const eckit::LocalConfiguration onlineDiag(config, "online diagnostics");
    bool writeinc = onlineDiag.getBool("write increment", false);

    if (writeinc) {
      // print log
      Log::info() << "Write Increment - starting: " << loop << std::endl << std::endl;

      const eckit::LocalConfiguration incConf(onlineDiag, "increment");

      // write increment
      dx.write(incConf);

      // print log
      Log::info() << std::endl << "Write Increment: done." << std::endl;
    }
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void writeKrylovBasis(const eckit::Configuration & config,
                        const ControlIncrement<MODEL, OBS> & dx,
                        const int & loop) {
// Write out the increment
  if (config.has("online diagnostics")) {
    eckit::LocalConfiguration diagConf(config, "online diagnostics");
    bool writeinc = diagConf.getBool("write basis", false);

    if (writeinc) {
      // print log
      Log::info() << "Write Krylov Basis: starting: " << loop << std::endl;

      eckit::LocalConfiguration basisConf(diagConf, "krylov basis");
      eckit::LocalConfiguration basisStateConf(basisConf, "state component");
      basisStateConf.set("iteration", loop);
      basisConf.set("state component", basisStateConf);

      // write increment
      dx.write(basisConf);

      // print log
      Log::info() << "Write Krylov Basis: done." << std::endl;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template<typename MODEL, typename OBS>
void writeEigenvectors(const eckit::Configuration & diagConf,
                       const std::vector<double> & diag,
                       const std::vector<double> & sub,
                       const std::vector<double> & rhs,
                       std::vector<std::unique_ptr<ControlIncrement<MODEL, OBS>>> & zvecs,
                       std::vector<std::unique_ptr<ControlIncrement<MODEL, OBS>>> & hvecs,
                       const HtRinvHMatrix<MODEL, OBS> & HtRinvH,
                       ControlIncrement<MODEL, OBS> & temp,
                       ControlIncrement<MODEL, OBS> & eigenv,
                       ControlIncrement<MODEL, OBS> & eigenz) {
  if (diagConf.has("online diagnostics")) {
    int maxEigen = diagConf.getInt("online diagnostics.max eigenvectors", 0);

    const double nn = rhs.size();
    Eigen::MatrixXd TT = Eigen::MatrixXd::Zero(nn, nn);

    for (int ii = 0; ii < nn; ++ii) {
      TT(ii, ii) = diag[ii];
      if (ii > 0) TT(ii-1, ii) = sub[ii-1];
      if (ii < nn - 1) TT(ii+1, ii) = sub[ii];
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> soluce(TT);
    Eigen::MatrixXd eigenvecT = soluce.eigenvectors().real();
    Eigen::VectorXd eigenvalT = soluce.eigenvalues().real();

    // Compute the eigenvectors (y = Zx)
    for (int ii = 0; ii < maxEigen && ii < nn; ++ii) {
      eigenz.zero();
      eigenv.zero();
      for (unsigned int jj = 0; jj < nn; ++jj) {
        temp.zero();
        temp = *zvecs[jj];
        temp *= eigenvecT.coeff(jj, nn - 1 - ii);
        eigenz += temp;
        temp.zero();
        temp = *hvecs[jj];
        temp *= eigenvecT.coeff(jj, nn - 1 - ii);
        eigenv += temp;
      }

      // Save the eigenvector
      if (diagConf.has("online diagnostics.eigenvector")) {
        eckit::LocalConfiguration basisConf(diagConf, "online diagnostics.eigenvector");
        eckit::LocalConfiguration basisStateConf(basisConf, "state component");
        basisStateConf.set("iteration", ii);
        basisConf.set("state component", basisStateConf);
        eigenz.write(basisConf);
      } else if (diagConf.has("online diagnostics.eigenvector to structured grid")) {
        eckit::LocalConfiguration eigenLatlonConf(diagConf,
                    "online diagnostics.eigenvector to structured grid");
        eigenLatlonConf.set("filename prefix",
              eigenLatlonConf.getString("filename prefix")+std::to_string(ii));
        // Eigenvector context has no meaningful State; therefore can't provide regular output on
        // pressure levels. Check here that model levels were requested:
        ASSERT(eigenLatlonConf.has("model levels") && !eigenLatlonConf.has("pressure levels"));
        const StructuredGridWriter<MODEL> latlon(eigenLatlonConf, eigenz.geometry());
        latlon.interpolateAndWrite(eigenz.state());
      }

      // Verification that eigenz is an eigenvector:
      // A.eigenv = eigenv + HtRinvH.Beigenv = eigenv + HtRinvH.eigenz = lambda eigenv
      temp.zero();
      HtRinvH.multiply(eigenz, temp);
      temp += eigenv;
      eigenv *= eigenvalT.coeff(nn - 1 - ii);
      temp -= eigenv;
      Log::info() << "Eigenvalue " << ii+1 << " : " << eigenvalT.coeff(nn - 1 - ii) << std::endl;
      Log::info() << "Norm A*y-lambda*y = " << dot_product(temp, temp) << std::endl;
      Log::info() << "Eigenvector " << ii+1 << " : " << eigenz << std::endl;

      Log::test() << "Eigenvalue " << ii+1 << " : " << eigenvalT.coeff(nn - 1 - ii) << std::endl;
      Log::test() << "Norm eigenvector = " << dot_product(eigenz, eigenz) << std::endl;
    }  // end for()
  }  // end if()
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_MINIMIZERUTILS_H_
