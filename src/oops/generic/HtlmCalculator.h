/*
 * (C) Copyright 2022-2023 UCAR.
 * (C) Crown copyright 2022-2023 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_HTLMCALCULATOR_H_
#define OOPS_GENERIC_HTLMCALCULATOR_H_

#include <string>
#include <vector>

#include "oops/base/IncrementSet.h"
#include "oops/generic/HtlmCalculatorCore.h"
#include "oops/generic/HtlmEnsemble.h"

namespace oops {

/*
 * Configuration options for HtlmCalculator:
 *
 * Keys:
 * ─────────────────────────────────────────────────────────────────────────────
 * "regularization"     : (Optional) configuration for regularization options
 *
 *  regularization sub keys
 *
 *    parts                : (Optional) specific regularization for specific lat lon
 *                                        locations, see HtlmRegularization.h for details
 *
 *    max condition number : (Optional) Maximum allowed condition number for SVD.
 *
 *    min singular value   : (Optional) Minimum singular value threshold.
 */

//------------------------------------------------------------------------------

/// \brief MODEL-templated front end for the HTLM coefficient calculation.
///
/// \details All of the dense linear algebra lives in HtlmCalculatorCore, which is not
/// templated and is compiled once. This class does the MODEL-dependent work: reading
/// sizes off the Geometry, obtaining the RMS values from the ensemble, and converting
/// IncrementSets into the atlas FieldSets the core consumes.
template <typename MODEL>
class HtlmCalculator {
  typedef Geometry<MODEL>                                          Geometry_;
  typedef HtlmEnsemble<MODEL>                                      HtlmEnsemble_;
  typedef Increment<MODEL>                                         Increment_;
  typedef IncrementSet<MODEL>                                      IncrementSet_;

 public:
  HtlmCalculator(const eckit::Configuration &,
                 const Variables &,
                 const Geometry_ &,
                 const atlas::idx_t,
                 const HtlmEnsemble_ &,
                 const std::vector<atlas::idx_t> &);
  static const std::string classname() {return "oops::HtlmCalculator";}
  void setOfCoeffs(const IncrementSet_ &, const IncrementSet_ &, atlas::FieldSet &) const;

 private:
  /// Component-dependent regularization needs a FieldSet shaped like the update variables;
  /// it is only built when the configuration asks for it. Returns an empty FieldSet otherwise.
  static atlas::FieldSet regularizationFieldSet(const eckit::Configuration &,
                                                const Geometry_ &,
                                                const Variables &);

  const Variables & updateVars_;
  const atlas::idx_t nLevels_;
  const atlas::idx_t ensembleSize_;
  HtlmCalculatorCore core_;
};

//------------------------------------------------------------------------------

template <typename MODEL>
atlas::FieldSet HtlmCalculator<MODEL>::regularizationFieldSet(
    const eckit::Configuration & config,
    const Geometry_ & updateGeometry,
    const Variables & updateVars) {
  atlas::FieldSet fset;
  if (config.getSubConfiguration("regularization").has("parts")) {
    Increment_ regularizationIncrement(updateGeometry, updateVars, util::DateTime());
    fset = regularizationIncrement.fieldSet().fieldSet();
  }
  return fset;
}

//------------------------------------------------------------------------------

template <typename MODEL>
HtlmCalculator<MODEL>::HtlmCalculator(const eckit::Configuration & config,
                                      const Variables & updateVars,
                                      const Geometry_ & updateGeometry,
                                      const atlas::idx_t influenceSize,
                                      const HtlmEnsemble_ & ensemble,
                                      const std::vector<atlas::idx_t> & owned)
: updateVars_(updateVars),
  nLevels_(updateGeometry.variableSizes(updateVars_)[0]),
  ensembleSize_(ensemble.size()),
  core_(config, updateVars_, nLevels_, influenceSize, ensembleSize_,
        ensemble.getRmsVals(updateVars_, nLevels_), owned,
        regularizationFieldSet(config, updateGeometry, updateVars_)) {}

//------------------------------------------------------------------------------

template<typename MODEL>
void HtlmCalculator<MODEL>::setOfCoeffs(const IncrementSet_ & linearEnsemble,
                                        const IncrementSet_ & linearErrors,
                                        atlas::FieldSet & coeffsFSet) const {
  std::vector<atlas::FieldSet> linearEnsembleFSets;
  std::vector<atlas::FieldSet> linearErrorsFSets;
  linearEnsembleFSets.reserve(ensembleSize_);
  linearErrorsFSets.reserve(ensembleSize_);
  for (atlas::idx_t m = 0; m < ensembleSize_; ++m) {
    linearEnsembleFSets.push_back(linearEnsemble[m].fieldSet().fieldSet());
    linearErrorsFSets.push_back(linearErrors[m].fieldSet().fieldSet());
  }
  core_.setOfCoeffs(linearEnsembleFSets, linearErrorsFSets, coeffsFSet);
}

//------------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_HTLMCALCULATOR_H_
