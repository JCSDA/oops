/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_RUNS_APPLICATION_H_
#define OOPS_RUNS_APPLICATION_H_

#include <iostream>
#include <string>

#include "eckit/mpi/Comm.h"

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

/// Parameters for testing in an oops application.
///
/// This is used only for further YAML validation, but is *not* used in setting up the tests in
/// oops::Run. See the more detailed warning inside `ApplicationParameters` below.
class ApplicationTestParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(ApplicationTestParameters, Parameters);

 public:
  OptionalParameter<std::string> mpiPattern{"mpi pattern", this};
  RequiredParameter<std::string> referenceFilename{"reference filename", this};

  OptionalParameter<double> floatAbsTol{"float absolute tolerance", this};
  OptionalParameter<double> floatRelTol{"float relative tolerance", this};
  OptionalParameter<int> integerTol{"integer tolerance", this};

  OptionalParameter<std::string> logOutputFilename{"log output filename", this};
  OptionalParameter<std::string> testOutputFilename{"test output filename", this};
};

// -----------------------------------------------------------------------------

/// Base class for top-level parameters of oops applications.
///
/// Parameters that all applications are expected to support should go here. Individual
/// applications should implement their own specialized parameters inherited from this base
/// class.
class ApplicationParameters : public Parameters {
  OOPS_ABSTRACT_PARAMETERS(ApplicationParameters, Parameters);

 public:
  /// Parameters used by regression tests comparing results produced by the application against
  /// known good outputs.
  ///
  /// \warning This is used only for further YAML validation, but is *not* used in setting up the
  /// tests in oops::Run. This is because oops::Run sets up the testing infrastructure *before* it
  /// executes the oops application, i.e., before we have an opportunity to validate the contents
  /// of the YAML via the Parameters. The reasons for adding a test Parameter to every application
  /// are,
  /// - to avoid duplicating code in the Parameters of every application
  /// - to work around the validation error the Parameters throws when it reads the test portion of
  ///   the YAML file if there is no test Parameter.
  /// - to add another layer of checks on the YAML's test portion, in case any errors were not
  ///   found by the parser in oops::Run.
  /// Long term, it would be more ideal to define the test Parameter at a level where it will
  /// directly control the reading of the test data.
  OptionalParameter<ApplicationTestParameters> test{"test", this};

  /// Output JSON Schema to a file specified by first argument.
  void outputSchema(const std::string & outputPath) const;
};

// -----------------------------------------------------------------------------

class Application : public util::Printable {
 public:
  explicit Application(const eckit::mpi::Comm & comm) : comm_(comm) {}
  virtual ~Application() {}
  virtual int execute(const eckit::Configuration &, bool validate) const = 0;

  /// This method aborts. Sub-class should override to output JSON schema.
  virtual void outputSchema(const std::string & outputPath) const;
  /// This method aborts. Sub-class should override to perform schema validation only
  virtual void validateConfig(const eckit::Configuration & fullConfig) const;

 protected:
  const eckit::mpi::Comm& getComm() const {return comm_;}

 private:
  const eckit::mpi::Comm & comm_;
  virtual std::string appname() const = 0;
  virtual void print(std::ostream & os) const {os << appname();}
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_RUNS_APPLICATION_H_
