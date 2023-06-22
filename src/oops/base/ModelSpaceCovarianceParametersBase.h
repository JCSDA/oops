/*
 * (C) Copyright 2020 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_MODELSPACECOVARIANCEPARAMETERSBASE_H_
#define OOPS_BASE_MODELSPACECOVARIANCEPARAMETERSBASE_H_

#include <string>
#include <vector>

#include "oops/interface/LinearVariableChange.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

namespace oops {

/// \brief Base class for classes storing parameters of a particular model-space error covariance
/// implementation.
template <typename MODEL>
class ModelSpaceCovarianceParametersBase : public Parameters {
  OOPS_ABSTRACT_PARAMETERS(ModelSpaceCovarianceParametersBase, Parameters)
 public:
  /// \brief Covariance model name.
  ///
  /// \note This parameter is marked as optional because it is only required in certain
  /// circumstances (e.g. when covariance model parameters are deserialized into a
  /// ModelSpaceCovarianceParametersWrapper and used by CovarianceFactory to instantiate
  /// a covariance model whose type is determined at runtime), but not others (e.g. in tests
  /// written with a particular model in mind). ModelSpaceCovarianceParametersWrapper will throw an
  /// exception if this parameter is not provided.
  typedef typename LinearVariableChange<MODEL>::Parameters_ Parameters_;
  OptionalParameter<std::string> covarianceModel{"covariance model", this};

  OptionalParameter<size_t> randomizationSize{"randomization size", this};
  Parameter<bool> fullInverse{"full inverse", false, this};
  Parameter<int> fullInverseIterations{"full inverse iterations", 10, this};
  Parameter<double> fullInverseAccuracy{"full inverse accuracy", 1.0e-3, this};
  OptionalParameter<Parameters_> variableChange{"linear variable change",
    this};
};

}  // namespace oops

#endif  // OOPS_BASE_MODELSPACECOVARIANCEPARAMETERSBASE_H_
