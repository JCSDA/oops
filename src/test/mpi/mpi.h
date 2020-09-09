/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_MPI_MPI_H_
#define TEST_MPI_MPI_H_

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "eckit/testing/Test.h"

#include "oops/../test/TestEnvironment.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/Expect.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace eckit
{
  // Don't use the contracted output for these types: the current implementation works only
  // with integer types.
  template <> struct VectorPrintSelector<util::DateTime> { typedef VectorPrintSimple selector; };
}  // namespace eckit

namespace test {

class TestParameters : public oops::Parameters {
 public:
  oops::RequiredParameter<std::vector<util::DateTime>> values{"values", this};
};

CASE("mpi/mpi/defaultCommunicators") {
  const eckit::mpi::Comm & world = oops::mpi::world();
  size_t worldsize = world.size();
  EXPECT_EQUAL(worldsize, 4);

  const eckit::mpi::Comm & talk_to_myself = oops::mpi::myself();
  size_t myownsize = talk_to_myself.size();
  EXPECT_EQUAL(myownsize, 1);
}

CASE("mpi/mpi/allGathervUsingSerialize") {
  const eckit::Configuration &conf = TestEnvironment::config();
  const eckit::mpi::Comm &comm = oops::mpi::world();

  TestParameters localParams;
  const size_t rank = comm.rank();
  localParams.deserialize(conf.getSubConfiguration("local" + std::to_string(rank)));
  const std::vector<util::DateTime> &localValues = localParams.values;

  TestParameters globalParams;
  globalParams.deserialize(conf.getSubConfiguration("global"));
  const std::vector<util::DateTime> &expectedGlobalValues = globalParams.values;

  size_t numGlobalValues;
  comm.allReduce(localValues.size(), numGlobalValues, eckit::mpi::Operation::SUM);

  std::vector<util::DateTime> globalValues(numGlobalValues);
  oops::mpi::allGathervUsingSerialize(comm, localValues.begin(), localValues.end(),
                                      globalValues.begin());
  EXPECT_EQUAL(globalValues, expectedGlobalValues);
}

class Mpi : public oops::Test {
 private:
  std::string testid() const override {return "test::mpi::mpi";}

  void register_tests() const override {}
};

}  // namespace test

#endif  // TEST_MPI_MPI_H_
