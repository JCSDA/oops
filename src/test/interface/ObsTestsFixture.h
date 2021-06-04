/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_INTERFACE_OBSTESTSFIXTURE_H_
#define TEST_INTERFACE_OBSTESTSFIXTURE_H_

#include <memory>
#include <vector>

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/ObsSpaces.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------
/// Fixture for observations-related tests
/// Gets created only once per test runs, ObsSpaces, configuration under `observations`
/// and window characteristics get saved
template <typename OBS>
class ObsTestsFixture : private boost::noncopyable {
  typedef oops::ObsSpaces<OBS>  ObsSpaces_;

 public:
  /// accessors to observation window
  static const util::DateTime & tbgn() {return *getInstance().tbgn_;}
  static const util::DateTime & tend() {return *getInstance().tend_;}
  /// accessor to a jj-th obs type config
  static eckit::LocalConfiguration & config(size_t jj) {return getInstance().configs_.at(jj);}
  /// accessor to a all obs spaces
  static ObsSpaces_ & obspace()        {return *getInstance().ospaces_;}
  static const eckit::mpi::Comm & comm()   {return getInstance().comm_;}

  static void reset() {
    obspace().save();
    getInstance().ospaces_.reset();
    getInstance().tend_.reset();
    getInstance().tbgn_.reset();
  }

 private:
  ObsTestsFixture(): comm_(oops::mpi::world()), tbgn_(), tend_(), ospaces_() {
    tbgn_.reset(new util::DateTime(TestEnvironment::config().getString("window begin")));
    tend_.reset(new util::DateTime(TestEnvironment::config().getString("window end")));
    configs_ = TestEnvironment::config().getSubConfigurations("observations");
    eckit::LocalConfiguration obsconfig =
           TestEnvironment::config().getSubConfiguration("observations");
    ospaces_.reset(new ObsSpaces_(obsconfig, comm_, *tbgn_, *tend_));
  }

  ~ObsTestsFixture() {}

  static ObsTestsFixture<OBS>& getInstance() {
    static ObsTestsFixture<OBS> theObsTestsFixture;
    return theObsTestsFixture;
  }

  const eckit::mpi::Comm & comm_;
  std::unique_ptr<const util::DateTime> tbgn_;
  std::unique_ptr<const util::DateTime> tend_;
  std::vector<eckit::LocalConfiguration> configs_;
  std::unique_ptr<ObsSpaces_> ospaces_;
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_OBSTESTSFIXTURE_H_
