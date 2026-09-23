/*
 * (C) Copyright 2026 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_HTLMCALCULATORCORE_H_
#define OOPS_GENERIC_HTLMCALCULATORCORE_H_

#include <memory>
#include <vector>

#include "atlas/field/FieldSet.h"
#include "atlas/library/config.h"

namespace eckit {
class Configuration;
}  // namespace eckit

namespace oops {

class Variables;

/// \brief Model-independent core of HtlmCalculator.
///
/// \details Owns the dense linear-algebra workspace (the preconditioned influence matrix,
/// its SVD, and the per-point linear-error vector) and performs the per-gridpoint singular
/// value decomposition and coefficient solve. None of that work depends on MODEL: it reads
/// and writes atlas FieldSets and operates on plain Eigen types.
class HtlmCalculatorCore {
 public:
  HtlmCalculatorCore(const eckit::Configuration & config,
                     const Variables & updateVars,
                     const atlas::idx_t nLevels,
                     const atlas::idx_t influenceSize,
                     const atlas::idx_t ensembleSize,
                     const atlas::FieldSet & rmsVals,
                     const std::vector<atlas::idx_t> & owned,
                     const atlas::FieldSet & regularizationFieldSet);
  ~HtlmCalculatorCore();

  HtlmCalculatorCore(const HtlmCalculatorCore &) = delete;
  HtlmCalculatorCore & operator=(const HtlmCalculatorCore &) = delete;

  // Compute the coefficients for every owned gridpoint and level, writing into coeffsFSet.
  // linearEnsemble and linearErrors are indexed by ensemble member.
  void setOfCoeffs(const std::vector<atlas::FieldSet> & linearEnsemble,
                   const std::vector<atlas::FieldSet> & linearErrors,
                   atlas::FieldSet & coeffsFSet) const;

 private:
  // The workspace is wrapped up in a subclass, so that it is not exposed by this header.
  // This keeps the very expensive Eigen::BDCSVD instantiation from being repeated for
  // each translation unit that includes HTLM code.
  class Workspace;
  std::unique_ptr<Workspace> impl_;
};

}  // namespace oops

#endif  // OOPS_GENERIC_HTLMCALCULATORCORE_H_
