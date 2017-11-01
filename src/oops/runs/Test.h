/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_RUNS_TEST_H_
#define OOPS_RUNS_TEST_H_

#include <cstring>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/tokenizer.hpp>

#include "util/Logger.h"
#include "oops/runs/Application.h"
#include "test/TestEnvironment.h"
#include "eckit/config/Configuration.h"

namespace oops {

// -----------------------------------------------------------------------------

class Test : public Application {
 public:
  Test() {}
  virtual ~Test() {}
  int execute(const eckit::Configuration & config) const;
 private:
  virtual void register_tests() const =0;
  virtual std::string testid() const =0;
  static bool init_unit_test() {return true;}
  std::string appname() const {return "oops::Test running " + testid();}
};

// -----------------------------------------------------------------------------

int Test::execute(const eckit::Configuration & config) const {
// Setup configuration for tests
  test::TestEnvironment::getInstance().setup(config);

// Extract the runtime config for the tests from the config file.
  std::string args = config.getString("test_framework_runtime_config");

// Create a string version of argv
  std::vector<std::string> argvec;
  argvec.push_back(std::string("abcd"));

  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> sep(" \n\t");
  tokenizer tok(args, sep);
  for (tokenizer::iterator it = tok.begin(); it != tok.end(); ++it) {
    argvec.push_back(*it);
  }

// Generate the argc and argv arguments for unit_test_main(...)
  int argc = argvec.size();
  char * argv[argc];
  for (int i = 0; i < argc; ++i) {
    argv[i] = new char[argvec[i].size()+1];
    strcpy(argv[i], argvec[i].c_str());
  }

// Run the tests
  Log::trace() << "Registering the unit tests" << std::endl;
  register_tests();
  Log::trace() << "Running the unit tests" << std::endl;
  int result = boost::unit_test::unit_test_main(&init_unit_test, argc, argv);
  Log::trace() << "Finished running the unit tests" << std::endl;
  Log::error() << "Finished running the unit tests, result = " << result << std::endl;

// Tidy up
  for (int i = 0; i < argc; ++i) {
    delete [] argv[i];
  }

// Return test status
  return result;
}

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_RUNS_TEST_H_
