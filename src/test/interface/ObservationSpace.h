/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_OBSERVATIONSPACE_H_
#define TEST_INTERFACE_OBSERVATIONSPACE_H_

#include <string>
#include <vector>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/runs/Test.h"
#include "oops/interface/ObservationSpace.h"
#include "test/TestEnvironment.h"
#include "eckit/config/LocalConfiguration.h"
#include "util/DateTime.h"

namespace test {

// -----------------------------------------------------------------------------
template <typename MODEL> class ObservationSpaceFixture : private boost::noncopyable {
  typedef oops::ObservationSpace<MODEL>        ObsSpace_;

 public:
  static const eckit::Configuration & config() {return *getInstance().conf_;}
  static const util::DateTime & bgn()  {return *getInstance().tbgn_;}
  static const util::DateTime & end()  {return *getInstance().tend_;}
  static const ObsSpace_ &  obspace()  {return *getInstance().obspace_;}

 private:
  static ObservationSpaceFixture<MODEL>& getInstance() {
    static ObservationSpaceFixture<MODEL> theObservationSpaceFixture;
    return theObservationSpaceFixture;
  }

  ObservationSpaceFixture() {
    tbgn_.reset(new util::DateTime(TestEnvironment::config().getString("window_begin")));
    tend_.reset(new util::DateTime(TestEnvironment::config().getString("window_end")));

    std::vector<eckit::LocalConfiguration> obsConfs;
    TestEnvironment::config().get("Observations", obsConfs);
    BOOST_CHECK(obsConfs.size() > 0);
    const eckit::LocalConfiguration obsConf(obsConfs[0], "Observation");
    obspace_.reset(new ObsSpace_(obsConf, *tbgn_, *tend_));
  }

  ~ObservationSpaceFixture() {}

  boost::scoped_ptr<const util::DateTime> tbgn_;
  boost::scoped_ptr<const util::DateTime> tend_;
  boost::scoped_ptr<ObsSpace_> obspace_;
};
// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef ObservationSpaceFixture<MODEL> Test_;

  BOOST_CHECK_EQUAL(Test_::obspace().windowStart(), Test_::bgn());
  BOOST_CHECK_EQUAL(Test_::obspace().windowEnd(),   Test_::end());
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testCopyConstructor() {
  typedef ObservationSpaceFixture<MODEL> Test_;
  typedef oops::ObservationSpace<MODEL>  ObsSpace_;

  boost::scoped_ptr<ObsSpace_> other(new ObsSpace_(Test_::obspace()));
  BOOST_CHECK(other.get());

  other.reset();
  BOOST_CHECK(!other.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> class ObservationSpace : public oops::Test {
 public:
  ObservationSpace() {}
  virtual ~ObservationSpace() {}
 private:
  std::string testid() const {return "test::ObservationSpace<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/ObservationSpace");

    ts->add(BOOST_TEST_CASE(&testConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testCopyConstructor<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSERVATIONSPACE_H_
