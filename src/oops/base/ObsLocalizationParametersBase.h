/*
 * (C) Copyright 2021- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

/// \brief Base class of classes storing parameters controlling specific
/// observation-space localizations.
class ObsLocalizationParametersBase : public Parameters {
  OOPS_ABSTRACT_PARAMETERS(ObsLocalizationParametersBase, Parameters)
 public:
  /// \brief Observation localization method.
  RequiredParameter<std::string> method{"localization method", this};
};

}  // namespace oops
