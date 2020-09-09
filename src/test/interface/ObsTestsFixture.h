/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_INTERFACE_OBSTESTSFIXTURE_H_
#define TEST_INTERFACE_OBSTESTSFIXTURE_H_

#include <memory>
#include <string>

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

template <typename OBS>
class ObsTestsFixture : private boost::noncopyable {
  typedef oops::ObsSpaces<OBS>  ObsSpaces_;

 public:
  static const util::DateTime & tbgn() {return *getInstance().tbgn_;}
  static const util::DateTime & tend() {return *getInstance().tend_;}
  static ObsSpaces_ & obspace()        {return *getInstance().ospaces_;}

 private:
  ObsTestsFixture(): tbgn_(), tend_(), ospaces_() {
    tbgn_.reset(new util::DateTime(TestEnvironment::config().getString("window begin")));
    tend_.reset(new util::DateTime(TestEnvironment::config().getString("window end")));
    ospaces_.reset(new ObsSpaces_(TestEnvironment::config(), oops::mpi::world(), *tbgn_, *tend_));
  }

  ~ObsTestsFixture() {}

  static ObsTestsFixture<OBS>& getInstance() {
    static ObsTestsFixture<OBS> theObsTestsFixture;
    return theObsTestsFixture;
  }

  std::unique_ptr<const util::DateTime> tbgn_;
  std::unique_ptr<const util::DateTime> tend_;
  std::unique_ptr<ObsSpaces_> ospaces_;
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_OBSTESTSFIXTURE_H_
