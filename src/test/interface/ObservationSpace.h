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

#include "eckit/config/LocalConfiguration.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/runs/Test.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  util::DateTime bgn(TestEnvironment::config().getString("window_begin"));
  util::DateTime end(TestEnvironment::config().getString("window_end"));

  const eckit::LocalConfiguration obsconf(TestEnvironment::config(), "Observations");
  std::vector<eckit::LocalConfiguration> conf;
  obsconf.get("ObsTypes", conf);

  for (std::size_t jj = 0; jj < conf.size(); ++jj) {
    oops::ObservationSpace<MODEL> obspace(conf[jj], bgn, end);

    BOOST_CHECK_EQUAL(obspace.windowStart(), bgn);
    BOOST_CHECK_EQUAL(obspace.windowEnd(),   end);
  }
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

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_OBSERVATIONSPACE_H_
