/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef TEST_INTERFACE_OBSOPERATOR_H_
#define TEST_INTERFACE_OBSOPERATOR_H_

#include <string>
#include <vector>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/runs/Test.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObservationSpace.h"
#include "test/TestEnvironment.h"
#include "eckit/config/LocalConfiguration.h"
#include "util/DateTime.h"

namespace test {

// -----------------------------------------------------------------------------
template <typename MODEL> class ObsOperFixture : private boost::noncopyable {
  typedef oops::ObservationSpace<MODEL>  ObsSpace_;

 public:
  static const ObsSpace_   & obspace() {return *getInstance().obspace_;}

 private:
  static ObsOperFixture<MODEL>& getInstance() {
    static ObsOperFixture<MODEL> theObsOperFixture;
    return theObsOperFixture;
  }

  ObsOperFixture() {
    const util::DateTime tbgn(TestEnvironment::config().getString("window_begin"));
    const util::DateTime tend(TestEnvironment::config().getString("window_end"));

    std::vector<eckit::LocalConfiguration> obsConfs;
    TestEnvironment::config().get("Observations", obsConfs);
    BOOST_CHECK(obsConfs.size() > 0);
    const eckit::LocalConfiguration obsConf(obsConfs[0], "Observation");
    obspace_.reset(new ObsSpace_(obsConf, tbgn, tend));
  }

  ~ObsOperFixture() {}

  boost::scoped_ptr<ObsSpace_>        obspace_;
};
// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef ObsOperFixture<MODEL> Test_;
  typedef oops::ObsOperator<MODEL>  ObsOperator_;

  boost::scoped_ptr<ObsOperator_> ov(new ObsOperator_(Test_::obspace()));
  BOOST_CHECK(ov.get());

  ov.reset();
  BOOST_CHECK(!ov.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> class ObsOperator : public oops::Test {
 public:
  ObsOperator() {}
  virtual ~ObsOperator() {}
 private:
  std::string testid() const {return "test::ObsOperator<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/ObsOperator");

    ts->add(BOOST_TEST_CASE(&testConstructor<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSOPERATOR_H_
