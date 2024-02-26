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
  static const util::TimeWindow & timeWindow() {return *getInstance().timeWindow_;}
  /// accessor to a jj-th obs type config
  static eckit::LocalConfiguration & config(size_t jj) {return getInstance().configs_.at(jj);}
  /// accessor to a all obs spaces
  static ObsSpaces_ & obspace()        {return *getInstance().ospaces_;}
  static const eckit::mpi::Comm & comm()   {return getInstance().comm_;}

  static void reset() {
    obspace().save();
    getInstance().ospaces_.reset();
    getInstance().timeWindow_.reset();
  }

 private:
  ObsTestsFixture(): comm_(oops::mpi::world()), timeWindow_(), ospaces_() {
    const eckit::LocalConfiguration conf(TestEnvironment::config());
    timeWindow_.reset(new util::TimeWindow(eckit::LocalConfiguration(conf, "time window")));
    configs_ = conf.getSubConfigurations("observations");
    eckit::LocalConfiguration obsconfig(conf, "observations");
    ospaces_.reset(new ObsSpaces_(obsconfig, comm_, *timeWindow_));
  }

  ~ObsTestsFixture() {}

  static ObsTestsFixture<OBS>& getInstance() {
    static ObsTestsFixture<OBS> theObsTestsFixture;
    return theObsTestsFixture;
  }

  const eckit::mpi::Comm & comm_;
  std::unique_ptr<const util::TimeWindow> timeWindow_;
  std::vector<eckit::LocalConfiguration> configs_;
  std::unique_ptr<ObsSpaces_> ospaces_;
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_OBSTESTSFIXTURE_H_
