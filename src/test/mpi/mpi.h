/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_MPI_MPI_H_
#define TEST_MPI_MPI_H_

#include <Eigen/Dense>

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

// -----------------------------------------------------------------------------------------------
class TestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(TestParameters, Parameters)
 public:
  oops::RequiredParameter<std::vector<util::DateTime>> datetime{"datetime", this};
  oops::RequiredParameter<std::vector<int>> int_{"int", this};
  oops::RequiredParameter<std::vector<std::string>> string{"string", this};
};

// -----------------------------------------------------------------------------------------------
template <typename T>
const std::vector<T> & getTestData(const TestParameters &params);

template <>
const std::vector<util::DateTime> & getTestData(const TestParameters &params) {
  return params.datetime;
}

template <>
const std::vector<int> & getTestData(const TestParameters &params) {
  return params.int_;
}

template <>
const std::vector<std::string> & getTestData(const TestParameters &params) {
  return params.string;
}

// -----------------------------------------------------------------------------------------------
template <typename T>
void testAllGatherv() {
  const eckit::Configuration &conf = TestEnvironment::config();
  const eckit::mpi::Comm &comm = oops::mpi::world();

  TestParameters localParams;
  const size_t rank = comm.rank();
  localParams.deserialize(conf.getSubConfiguration("local" + std::to_string(rank)));
  std::vector<T> values = getTestData<T>(localParams);

  TestParameters globalParams;
  globalParams.deserialize(conf.getSubConfiguration("global"));
  const std::vector<T> &expectedResult = getTestData<T>(globalParams);

  oops::mpi::allGatherv(comm, values);
  EXPECT_EQUAL(values, expectedResult);
}

// -----------------------------------------------------------------------------------------------
CASE("mpi/mpi/defaultCommunicators") {
  const eckit::mpi::Comm & world = oops::mpi::world();
  size_t worldsize = world.size();
  EXPECT_EQUAL(worldsize, 4);

  const eckit::mpi::Comm & talk_to_myself = oops::mpi::myself();
  size_t myownsize = talk_to_myself.size();
  EXPECT_EQUAL(myownsize, 1);
}
// -----------------------------------------------------------------------------------------------
CASE("mpi/mpi/allGathervUsingSerialize") {
  const eckit::Configuration &conf = TestEnvironment::config();
  const eckit::mpi::Comm &comm = oops::mpi::world();

  TestParameters localParams;
  const size_t rank = comm.rank();
  localParams.deserialize(conf.getSubConfiguration("local" + std::to_string(rank)));
  const std::vector<util::DateTime> &localValues = localParams.datetime;

  TestParameters globalParams;
  globalParams.deserialize(conf.getSubConfiguration("global"));
  const std::vector<util::DateTime> &expectedGlobalValues = globalParams.datetime;

  size_t numGlobalValues;
  comm.allReduce(localValues.size(), numGlobalValues, eckit::mpi::Operation::SUM);

  std::vector<util::DateTime> globalValues(numGlobalValues);
  oops::mpi::allGathervUsingSerialize(comm, localValues.begin(), localValues.end(),
                                      globalValues.begin());
  EXPECT_EQUAL(globalValues, expectedGlobalValues);
}
// -----------------------------------------------------------------------------------------------
CASE("mpi/mpi/allGathervInt") {
  testAllGatherv<int>();
}
// -----------------------------------------------------------------------------------------------
CASE("mpi/mpi/allGathervDateTime") {
  testAllGatherv<util::DateTime>();
}
// -----------------------------------------------------------------------------------------------
CASE("mpi/mpi/allGathervInt") {
  testAllGatherv<std::string>();
}
// -----------------------------------------------------------------------------------------------
CASE("mpi/mpi/SendReceive") {
  const eckit::Configuration &conf = TestEnvironment::config();
  const eckit::mpi::Comm &comm = oops::mpi::world();
  const size_t rank = comm.rank();
  int source = (rank + 3) % comm.size();
  int destination = (rank + 1) % comm.size();
  int tag_send = destination;
  int tag_recv = rank;

  util::DateTime sendValue(conf.getString("send"+ std::to_string(rank)));
  util::DateTime expectedValue(conf.getString("expected"+ std::to_string(rank)));
  util::DateTime receivedValue;

  if (rank < 3) {
    oops::mpi::send(comm, sendValue, destination, tag_send);
    oops::mpi::receive(comm, receivedValue, source, tag_recv);
  } else {
    oops::mpi::receive(comm, receivedValue, source, tag_recv);
    oops::mpi::send(comm, sendValue, destination, tag_send);
  }
  EXPECT_EQUAL(receivedValue, expectedValue);
}
// -----------------------------------------------------------------------------------------------
CASE("mpi/mpi/gatherSerializable") {
  const eckit::Configuration &conf = TestEnvironment::config();
  const eckit::mpi::Comm &comm = oops::mpi::world();

  TestParameters localParams;
  const size_t rank = comm.rank();
  localParams.deserialize(conf.getSubConfiguration("local" + std::to_string(rank)));
  const std::vector<util::DateTime> &localValues = localParams.datetime;

  TestParameters globalParams;
  globalParams.deserialize(conf.getSubConfiguration("global"));
  const std::vector<util::DateTime> &expectedGlobalValues = globalParams.datetime;

  size_t numGlobalValues;
  comm.allReduce(localValues.size(), numGlobalValues, eckit::mpi::Operation::SUM);

  std::vector<util::DateTime> globalValues(numGlobalValues);

  util::DateTime zeroDate("0001-01-01T00:00:00Z");
  for (size_t ii = 0; ii < numGlobalValues; ++ii) {
    globalValues[ii] = zeroDate;
  }

  std::vector<util::DateTime> zeroValues = globalValues;

  size_t root_gather = conf.getInt("root for gathering", 0);

  oops::mpi::gather(comm, localValues, globalValues, root_gather);
  if (rank == root_gather) {
    EXPECT_EQUAL(globalValues, expectedGlobalValues);
  } else {
    EXPECT_EQUAL(globalValues, zeroValues);
  }
}
// -----------------------------------------------------------------------------------------------
CASE("mpi/mpi/gatherDouble") {
  const eckit::Configuration &conf = TestEnvironment::config();
  const eckit::mpi::Comm &comm = oops::mpi::world();
  const size_t rank = comm.rank();

  std::vector<double> localDouble;
  conf.get("localDouble" + std::to_string(rank), localDouble);

  std::vector<double> globalDoubleExpected;
  conf.get("globalDouble", globalDoubleExpected);

  size_t numGlobalDouble;
  comm.allReduce(localDouble.size(), numGlobalDouble, eckit::mpi::Operation::SUM);

  std::vector<double> globalDouble(numGlobalDouble, 0.0);
  std::vector<double> zerosDouble = globalDouble;

  size_t root_gather = conf.getInt("root for gathering", 0);

  oops::mpi::gather(comm, localDouble, globalDouble, root_gather);

  if (rank == root_gather) {
    EXPECT_EQUAL(globalDouble, globalDoubleExpected);
  } else {
    EXPECT_EQUAL(globalDouble, zerosDouble);
  }
}
// -----------------------------------------------------------------------------------------------
CASE("mpi/mpi/allGatherEigen") {
  const eckit::mpi::Comm &comm = oops::mpi::world();
  const size_t rank = comm.rank();

  Eigen::VectorXd localEigen = rank * Eigen::VectorXd::Ones(4);

  Eigen::MatrixXd globalEigen(4, 4);
  globalEigen << Eigen::VectorXd::Zero(4),
                 Eigen::VectorXd::Zero(4),
                 Eigen::VectorXd::Zero(4),
                 Eigen::VectorXd::Zero(4);

  Eigen::MatrixXd expectedEigen(4, 4);
  expectedEigen << 0*Eigen::VectorXd::Ones(4),
                   1*Eigen::VectorXd::Ones(4),
                   2*Eigen::VectorXd::Ones(4),
                   3*Eigen::VectorXd::Ones(4);

  oops::mpi::allGather(comm, localEigen, globalEigen);
  EXPECT_EQUAL(expectedEigen, globalEigen);
}
// -----------------------------------------------------------------------------------------------
CASE("mpi/mpi/exclusiveScan") {
  const eckit::mpi::Comm &comm = oops::mpi::world();
  const size_t rank = comm.rank();

  size_t expectedResult = 0;
  for (size_t lowerRank = 0; lowerRank < rank; ++lowerRank)
    expectedResult += lowerRank;

  size_t result = rank;
  oops::mpi::exclusiveScan(comm, result);
  EXPECT_EQUAL(result, expectedResult);
}
// -----------------------------------------------------------------------------------------------
CASE("mpi/mpi/broadcastBool") {
  const eckit::mpi::Comm &comm = oops::mpi::world();
  const size_t rank = comm.rank();

  bool expectedResult = true;
  bool result;
  if (rank == 0) {
    result = expectedResult;
  } else {
    result = false;
  }
  oops::mpi::broadcastBool(comm, result, 0);
  EXPECT_EQUAL(result, expectedResult);
}
// -----------------------------------------------------------------------------------------------
CASE("mpi/mpi/broadcastString") {
  const eckit::mpi::Comm &comm = oops::mpi::world();
  const size_t rank = comm.rank();

  std::string expectedResult("Hello World");
  std::string result;
  if (rank == 0) {
    result = expectedResult;
  } else {
    result = "NULL";
  }
  oops::mpi::broadcastString(comm, result, 0);
  EXPECT_EQUAL(result, expectedResult);
}
// -----------------------------------------------------------------------------------------------
CASE("mpi/mpi/broadcastVector<int>") {
  const eckit::mpi::Comm &comm = oops::mpi::world();
  const size_t rank = comm.rank();

  std::vector<int> expectedResult = { 0, 1, 2, 3, 4, 5 };
  std::vector<int> result;
  if (rank == 0) {
    result = expectedResult;
  } else {
    result = { 10, 20 , 30 , 40 };
  }
  oops::mpi::broadcastVector<int>(comm, result, 0);
  EXPECT_EQUAL(result, expectedResult);
}
// -----------------------------------------------------------------------------------------------
CASE("mpi/mpi/broadcastVector<char>") {
  const eckit::mpi::Comm &comm = oops::mpi::world();
  const size_t rank = comm.rank();

  std::vector<char> expectedResult = { 'a' , 'b', 'c' };
  std::vector<char> result;
  if (rank == 0) {
    result = expectedResult;
  } else {
    result = { 'x', 'y' , 'z' };
  }
  oops::mpi::broadcastVector<char>(comm, result, 0);
  EXPECT_EQUAL(result, expectedResult);
}
// -----------------------------------------------------------------------------------------------
CASE("mpi/mpi/sendReceiveString") {
  const eckit::mpi::Comm &comm = oops::mpi::world();
  const size_t rank = comm.rank();

  std::string expectedResult("Hello World");
  std::string result;
  if (rank == 0) {
    result = expectedResult;
    for (size_t i = 1; i < comm.size(); ++i) {
        oops::mpi::sendString(comm, result, i);
    }
  } else {
    result = "NULL";
    oops::mpi::receiveString(comm, result, 0);
  }
  EXPECT_EQUAL(result, expectedResult);
}
// -----------------------------------------------------------------------------------------------

class Mpi : public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "test::mpi::mpi";}

  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test

#endif  // TEST_MPI_MPI_H_
