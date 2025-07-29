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

namespace test {
  class DummyApp: public oops::Application {
   public:
    explicit DummyApp(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}
    virtual ~DummyApp() = default;
    int execute(const eckit::Configuration & fullConfig) const override {
      std::string hello_str = fullConfig.getString("hello");
      oops::Log::info() << "hello " << hello_str << std::endl;
      return 0;
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
