/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/generic/PseudoModel.h"

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

namespace oops {

// -----------------------------------------------------------------------------

MultipleStateTemplateConfigValidator::MultipleStateTemplateConfigValidator(
    eckit::LocalConfiguration stateTemplateConf) {
  // template to define a generic state
  config.set("template", stateTemplateConf.getSubConfiguration("template"));
  // Starting state datetime for the first state
  config.set("start datetime", stateTemplateConf.getString("start datetime"));
  // total number of states in the list after excluding certain values
  config.set("number of states", stateTemplateConf.getUnsigned("number of states"));

  // pattern to be replaced by a padded integer for each state
  if (stateTemplateConf.has("pattern")) {
    config.set("pattern", stateTemplateConf.getString("pattern"));
  }

  // starting state padded integer value for the first state
  config.set("start", size_t(1));
  if (stateTemplateConf.has("start")) {
    config.set("start", stateTemplateConf.getUnsigned("start"));
  }
  // list of excluded state pattern values
  config.set("except", std::vector<size_t>{});
  if (stateTemplateConf.has("except")) {
    config.set("except", stateTemplateConf.getUnsignedVector("except"));
  }

  // zero padding
  config.set("zero padding", size_t(0));
  if (stateTemplateConf.has("zero padding")) {
    config.set("zero padding", stateTemplateConf.getUnsigned("zero padding"));
  }
}

// -----------------------------------------------------------------------------

eckit::LocalConfiguration MultipleStateTemplateConfigValidator::indexToStateConfig(
    const size_t & stateConfIndex, const util::Duration tstep) const {
  // Get correct index
  size_t patternValue = this->config.getUnsigned("start");
  util::DateTime dateTimeValue(this->config.getString("start datetime"));
  std::vector<size_t> except(this->config.getUnsignedVector("except"));
  for (size_t jj = 0; jj <= stateConfIndex; ++jj) {
    // Check for excluded states
    while (std::count(except.begin(), except.end(), patternValue)) {
      patternValue += 1;
    }
    // Update patternValue
    if (jj < stateConfIndex) {
      patternValue += 1;
      dateTimeValue += tstep;
    }
  }

  // Copy and update template configuration with the current pattern value and date
  eckit::LocalConfiguration stateConf(this->config.getSubConfiguration("template"));
  stateConf.set("date", dateTimeValue.toString());
  if (this->config.has("pattern")) {
    std::string pattern(this->config.getString("pattern"));
    util::seekAndReplace(stateConf, pattern, patternValue,
                         this->config.getUnsigned("zero padding"));
  }

  return stateConf;
}

// -----------------------------------------------------------------------------

MultipleStateConfigValidator::MultipleStateConfigValidator(
    eckit::LocalConfiguration multipleStateConf) {
  // list of the states in the grouping
  if (multipleStateConf.has("states")) {
    config.set("states", multipleStateConf.getSubConfigurations("states"));
  }
  // template to define states of the list of multiple states
  if (multipleStateConf.has("states from template")) {
    config.set("states from template",
        multipleStateConf.getSubConfiguration("states from template"));
  }
  // run identifier
  config.set("ID", 0);
  if (multipleStateConf.has("ID")) {
    config.set("ID", multipleStateConf.getUnsigned("ID"));
  }

  if (!config.has("states") && !config.has("states from template")) {
    throw eckit::UserError(
        "MultipleStateConfigValidator: both 'states' and 'states from template' are missing",
        Here());
  }
  if (config.has("states") && config.has("states from template")) {
    throw eckit::UserError(
        "MultipleStateConfigValidator: both 'states' and 'states from template' are present",
        Here());
  }
}

// -----------------------------------------------------------------------------

size_t MultipleStateConfigValidator::size() const {
  if (this->config.has("states")) {
    return this->config.getSubConfigurations("states").size();
  }
  auto templateConfig = this->config.getSubConfiguration("states from template");
  return templateConfig.getUnsigned("number of states");
}

// -----------------------------------------------------------------------------

eckit::LocalConfiguration MultipleStateConfigValidator::getStateConfig(
    const size_t & stateConfIndex, const util::Duration tstep) const {
  // Check number of states
  if (stateConfIndex >= this->size()) {
    throw eckit::UserError(
        "MultipleStateConfigValidator::getStateConfig index argument is too large",
        Here());
  }

  if (this->config.has("states")) {
    // Explicit states
    auto statesConfigs = this->config.getSubConfigurations("states");
    return statesConfigs[stateConfIndex];
  }
  MultipleStateTemplateConfigValidator validator{
    this->config.getSubConfiguration("states from template")
  };
  return validator.indexToStateConfig(stateConfIndex, tstep);
}

// -----------------------------------------------------------------------------

PseudoModelConfigValidator::PseudoModelConfigValidator(
    eckit::LocalConfiguration pseudoModelConf) {
  // list of the states in the grouping
  if (pseudoModelConf.has("states")) {
    config.set("states", pseudoModelConf.getSubConfigurations("states"));
  }
  // template to define states of the list of multiple states
  if (pseudoModelConf.has("states from template")) {
    config.set("states from template",
      pseudoModelConf.getSubConfiguration("states from template"));
  }
  // list of the runs in the Pseudo-model
  if (pseudoModelConf.has("multiple runs")) {
    config.set("multiple runs",
      pseudoModelConf.getSubConfigurations("multiple runs"));
  }

  // Name of the model type
  if (pseudoModelConf.getString("name") != "PseudoModel") {
      throw eckit::UserError(
          "PseudoModelConfigValidator: 'name' is not 'PseudoModel' for this PseudoModel"
          " configuration.",
          Here());
  }
  config.set("name", "PseudoModel");
  // time step between the states in a run
  config.set("tstep", pseudoModelConf.getString("tstep"));

  if (!config.has("states")
      && !config.has("states from template")
      && !config.has("multiple runs")) {
    throw eckit::UserError(
        "PseudoModelConfigValidator: 'states', 'states from template', and 'multiple runs'"
        " parameters are all missing. At least one must be specified.",
        Here());
  }
  if (config.has("states") && config.has("states from template")) {
    throw eckit::UserError(
        "PseudoModelConfigValidator: both states and states from template parameters are"
        " present. Only one or the other should be specified.",
        Here());
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops
