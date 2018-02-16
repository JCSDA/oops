/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef TEST_INTERFACE_OBSTESTSFIXTURE_H_
#define TEST_INTERFACE_OBSTESTSFIXTURE_H_

#include <string>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/runs/Test.h"
#include "oops/base/ObsSpaces.h"
#include "test/TestEnvironment.h"
#include "eckit/config/LocalConfiguration.h"
#include "util/DateTime.h"
#include "util/Logger.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsTestsFixture : private boost::noncopyable {
  typedef oops::ObsSpaces<MODEL>  ObsSpaces_;

 public:
  static const util::DateTime & tbgn() {return *getInstance().tbgn_;}
  static const util::DateTime & tend() {return *getInstance().tend_;}
  static ObsSpaces_ & obspace()        {return *getInstance().ospaces_;}

 private:
  ObsTestsFixture(): tbgn_(), tend_(), ospaces_() {
    tbgn_.reset(new util::DateTime(TestEnvironment::config().getString("window_begin")));
    tend_.reset(new util::DateTime(TestEnvironment::config().getString("window_end")));

    const eckit::LocalConfiguration conf(TestEnvironment::config(), "Observations");
    ospaces_.reset(new ObsSpaces_(conf, *tbgn_, *tend_));
  }

  ~ObsTestsFixture() {}

  static ObsTestsFixture<MODEL>& getInstance() {
    static ObsTestsFixture<MODEL> theObsTestsFixture;
    return theObsTestsFixture;
  }

  boost::scoped_ptr<const util::DateTime> tbgn_;
  boost::scoped_ptr<const util::DateTime> tend_;
  boost::scoped_ptr<ObsSpaces_> ospaces_;
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_OBSTESTSFIXTURE_H_
