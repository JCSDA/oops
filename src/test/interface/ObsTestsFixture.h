/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_INTERFACE_OBSTESTSFIXTURE_H_
#define TEST_INTERFACE_OBSTESTSFIXTURE_H_

#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/ObsSpaces.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
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
  static const eckit::mpi::Comm & comm()   {return *getCommPointerInstance();}

  static void reset() {
    obspace().save();
    getInstance().ospaces_.reset();
    getInstance().timeWindow_.reset();
  }

  /// \brief Set the communicator to be used by all obs spaces.
  ///
  /// To have an effect, this function must be called before the first call to any other member
  /// function except comm().
  static void setComm(const eckit::mpi::Comm &comm) {
    getCommPointerInstance() = &comm;
  }

 private:
  ObsTestsFixture(): timeWindow_(), ospaces_() {
    const eckit::LocalConfiguration conf(TestEnvironment::config());
    timeWindow_.reset(new util::TimeWindow(eckit::LocalConfiguration(conf, "time window")));

    // Many test YAMLs write "observations:" as a bare sequence of observers (sequence form),
    // whereas ObsSpaces expects the "observations" section, which may also carry settings
    // applying to every obs space (map form). Accept both, normalising the sequence form
    // into the mapping form so that those settings have somewhere to live.
    const bool haveObservers = conf.has("observations.observers");
    configs_ = haveObservers ? conf.getSubConfigurations("observations.observers")
                             : conf.getSubConfigurations("observations");

    eckit::LocalConfiguration obsconfig;
    if (haveObservers) {
      obsconfig = eckit::LocalConfiguration(conf, "observations");
    } else {
      oops::Log::info() << "WARNING: ObsTestsFixture: YAML 'observations' section with a "
              << "list of individual 'obs space' specs is a deprecated format. " << std::endl
              << "WARNING: Please nest the list of 'obs space' specs under an "
              << "'observations.observers:' key " << std::endl;
      obsconfig.set("observers", configs_);
      // Transform the sequence form into the mapping form for the subsequent
      // ObsSpaces constructor. Temporarily allow for the "obs data container"
      // option to be set at the top level of the test YAML (ie, sibling to
      // "time window" and "observations"), and used in the new map form.
      // TODO(srh): remove this once the test YAMLs have moved to the map form.
      std::string obsDataContainer;
      if (conf.get("obs data container", obsDataContainer))
        obsconfig.set("obs data container", obsDataContainer);
    }

    ospaces_.reset(new ObsSpaces_(obsconfig, *getCommPointerInstance(), *timeWindow_));
  }

  ~ObsTestsFixture() {}

  static ObsTestsFixture<OBS>& getInstance() {
    static ObsTestsFixture<OBS> theObsTestsFixture;
    return theObsTestsFixture;
  }

  static const eckit::mpi::Comm *& getCommPointerInstance() {
    static const eckit::mpi::Comm * theCommPointer = &oops::mpi::world();
    return theCommPointer;
  }

  std::unique_ptr<const util::TimeWindow> timeWindow_;
  std::vector<eckit::LocalConfiguration> configs_;
  std::unique_ptr<ObsSpaces_> ospaces_;
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_OBSTESTSFIXTURE_H_
