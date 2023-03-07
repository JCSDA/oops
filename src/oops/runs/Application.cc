/*
 * Copyright (C) British Crown (Met Office) & Contributors.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/runs/Application.h"

#include <fstream>

#include "oops/util/abor1_cpp.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Output JSON Schema to a file specified by the outputPath argument.
void ApplicationParameters::outputSchema(const std::string & outputPath) const {
  std::ofstream stream(outputPath);
  stream << jsonSchema().toString(true);
  stream.close();
}

// -----------------------------------------------------------------------------

/// This method aborts. Sub-class should override to output JSON schema.
void Application::outputSchema(const std::string & outputPath) const {
  ABORT(appname() + "::outputSchema not implemented");
}

/// This method aborts. Sub-class should override to perform schema validation only.
void Application::validateConfig(const eckit::Configuration & fullConfig) const {
  ABORT(appname() + "::validateConfig not implemented");
}

// -----------------------------------------------------------------------------

}  // namespace oops
