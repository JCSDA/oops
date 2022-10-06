/*
 * (C) Copyright 2022 British Crown (Met Office) & Contributors.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 *
 */
#include <string>

#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/runs/Run.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace test {
  class DummyAppParameters : public oops::ApplicationParameters {
    OOPS_CONCRETE_PARAMETERS(DummyAppParameters, ApplicationParameters);
   public:
    oops::RequiredParameter<std::string> hello{"hello", "Who to greet?", this};
  };
  class DummyApp: public oops::Application {
   public:
    explicit DummyApp(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}
    virtual ~DummyApp() = default;
    int execute(const eckit::Configuration & fullConfig, bool validate) const override {
      DummyAppParameters params;
      if (validate) params.validate(fullConfig);
      params.deserialize(fullConfig);
      std::string hello_str = params.hello;
      oops::Log::info() << "hello " << hello_str << std::endl;
      return 0;
    }
    void outputSchema(const std::string & outputPath) const override {
      DummyAppParameters params;
      params.outputSchema(outputPath);
    }
    void validateConfig(const eckit::Configuration & fullConfig) const override {
      DummyAppParameters params;
      params.validate(fullConfig);
    }
   private:
    std::string appname() const override {
      return "test::DummyApp";
    }
  };
}  // namespace test

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::DummyApp dummyApp;
  return run.execute(dummyApp);
}
