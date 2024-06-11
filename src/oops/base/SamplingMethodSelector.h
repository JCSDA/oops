/*
 * (C) Crown copyright 2023, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_SAMPLINGMETHODSELECTOR_H_
#define OOPS_BASE_SAMPLINGMETHODSELECTOR_H_

namespace oops {

class Variable;

/// \brief Maps each model variable and ObsDiagnostic to the observation location sampling method
/// producing the set of paths along which this variable should be interpolated.
class SamplingMethodSelector {
 public:
  virtual ~SamplingMethodSelector() {}

  /// \brief Returns the index of the sampling method to be used for the variable `varName`.
  virtual size_t methodIndex(const Variable &varName) const = 0;
};

}  // namespace oops

#endif  // OOPS_BASE_SAMPLINGMETHODSELECTOR_H_
