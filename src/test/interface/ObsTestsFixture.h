/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef TEST_INTERFACE_OBSTESTSFIXTURE_H_
#define TEST_INTERFACE_OBSTESTSFIXTURE_H_

#include <string>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/ObsOperators.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsTestsFixture : private boost::noncopyable {
  typedef oops::ObsOperators<MODEL>  ObsOperators_;

 public:
  static const util::DateTime & tbgn() {return *getInstance().tbgn_;}
  static const util::DateTime & tend() {return *getInstance().tend_;}
  static ObsOperators_ & hoper()       {return *getInstance().hoper_;}

 private:
  ObsTestsFixture(): tbgn_(), tend_(), hoper_() {
    tbgn_.reset(new util::DateTime(TestEnvironment::config().getString("window_begin")));
    tend_.reset(new util::DateTime(TestEnvironment::config().getString("window_end")));

    const eckit::LocalConfiguration conf(TestEnvironment::config(), "Observations");
    hoper_.reset(new ObsOperators_(conf, *tbgn_, *tend_));
  }

  ~ObsTestsFixture() {}

  static ObsTestsFixture<MODEL>& getInstance() {
    static ObsTestsFixture<MODEL> theObsTestsFixture;
    return theObsTestsFixture;
  }

  boost::scoped_ptr<const util::DateTime> tbgn_;
  boost::scoped_ptr<const util::DateTime> tend_;
  boost::scoped_ptr<ObsOperators_> hoper_;
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_OBSTESTSFIXTURE_H_
