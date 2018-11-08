/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_OBSVECTOR_H_
#define TEST_INTERFACE_OBSVECTOR_H_

#include <string>
#include <vector>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/shared_ptr.hpp>

#include "oops/interface/ObservationSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/runs/Test.h"
#include "oops/util/dot_product.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsVecFixture : private boost::noncopyable {
  typedef oops::ObservationSpace<MODEL>  ObsSpace_;

 public:
  static std::vector<boost::shared_ptr<ObsSpace_> > & obspace() {return getInstance().ospaces_;}

 private:
  static ObsVecFixture<MODEL>& getInstance() {
    static ObsVecFixture<MODEL> theObsVecFixture;
    return theObsVecFixture;
  }

  ObsVecFixture(): ospaces_() {
    util::DateTime bgn((TestEnvironment::config().getString("window_begin")));
    util::DateTime end((TestEnvironment::config().getString("window_end")));

    const eckit::LocalConfiguration obsconf(TestEnvironment::config(), "Observations");
    std::vector<eckit::LocalConfiguration> conf;
    obsconf.get("ObsTypes", conf);

    for (std::size_t jj = 0; jj < conf.size(); ++jj) {
      boost::shared_ptr<ObsSpace_> tmp(new ObsSpace_(conf[jj], bgn, end));
      ospaces_.push_back(tmp);
    }
  }

  ~ObsVecFixture() {}

  std::vector<boost::shared_ptr<ObsSpace_> > ospaces_;
};

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef ObsVecFixture<MODEL>  Test_;
  typedef oops::ObsVector<MODEL>  ObsVector_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    boost::scoped_ptr<ObsVector_> ov(new ObsVector_(*Test_::obspace()[jj]));
    BOOST_CHECK(ov.get());

    ov.reset();
    BOOST_CHECK(!ov.get());
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testCopyConstructor() {
  typedef ObsVecFixture<MODEL>  Test_;
  typedef oops::ObsVector<MODEL>  ObsVector_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    boost::scoped_ptr<ObsVector_> ov(new ObsVector_(*Test_::obspace()[jj]));

    boost::scoped_ptr<ObsVector_> other(new ObsVector_(*ov));
    BOOST_CHECK(other.get());

    other.reset();
    BOOST_CHECK(!other.get());

    BOOST_CHECK(ov.get());
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testNotZero() {
  typedef ObsVecFixture<MODEL>  Test_;
  typedef oops::ObsVector<MODEL>  ObsVector_;
  const double zero = 0.0;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ObsVector_ ov(*Test_::obspace()[jj]);

    ov.random();

    const double ovov2 = dot_product(ov, ov);
    BOOST_CHECK(ovov2 > zero);

    ov.zero();

    const double zz = dot_product(ov, ov);
    BOOST_CHECK(zz == zero);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> class ObsVector : public oops::Test {
 public:
  ObsVector() {}
  virtual ~ObsVector() {}
 private:
  std::string testid() const {return "test::ObsVector<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/ObsVector");

    ts->add(BOOST_TEST_CASE(&testConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testCopyConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testNotZero<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSVECTOR_H_
