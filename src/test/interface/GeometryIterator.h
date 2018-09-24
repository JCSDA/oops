/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_GEOMETRYITERATOR_H_
#define TEST_INTERFACE_GEOMETRYITERATOR_H_

#include <math.h>
#include <string>
#include <vector>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/GridPoint.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/interface/State.h"
#include "oops/runs/Test.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------
template <typename MODEL> class GeometryIteratorFixture : private boost::noncopyable {
 public:
  static const eckit::Configuration & getConfig() {return *getInstance().conf_;}

 private:
  static GeometryIteratorFixture<MODEL>& getInstance() {
    static GeometryIteratorFixture<MODEL> theGeometryIteratorFixture;
    return theGeometryIteratorFixture;
  }

  GeometryIteratorFixture() {
    conf_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "GeometryIterator"));
  }

  ~GeometryIteratorFixture() {}

  boost::scoped_ptr<const eckit::LocalConfiguration> conf_;
};
// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef oops::GeometryIterator<MODEL>        GeometryIterator_;
  typedef oops::Geometry<MODEL>        Geometry_;

  const eckit::LocalConfiguration
       geomConfig(GeometryIteratorFixture<MODEL>::getConfig(), "Geometry");
  Geometry_ geom(geomConfig);

  boost::scoped_ptr<GeometryIterator_> iter(new GeometryIterator_(geom.begin()));
  BOOST_CHECK(iter.get());

  iter.reset();
  BOOST_CHECK(!iter.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testIterator() {
  typedef oops::Geometry<MODEL>          Geometry_;
  typedef oops::GeometryIterator<MODEL>  GeometryIterator_;
  typedef oops::State<MODEL>             State_;

  const eckit::LocalConfiguration
       geomConfig(GeometryIteratorFixture<MODEL>::getConfig(), "Geometry");
  const eckit::LocalConfiguration stateConfig(GeometryIteratorFixture<MODEL>::getConfig(), "State");

  const oops::Variables vars(stateConfig);

  Geometry_ geom(geomConfig);
  State_ state(geom, vars, stateConfig);

  const double rms_conf = stateConfig.getDouble("rms");
  const double tol = stateConfig.getDouble("tolerance");

  double rms = 0;
  int n = 0;
  for (GeometryIterator_ i = geom.begin(); i != geom.end(); ++i, ++n) {
    oops::GridPoint gp = state.getPoint(i);
    oops::Log::debug() << *i << gp << std::endl;
    std::vector<double> vals = gp.getVals();
    rms += std::inner_product(vals.begin(), vals.end(), vals.begin(), 0.) /
                  vals.size();
  }
  rms = sqrt(rms/n);

  oops::Log::debug() << n << " rms from iterator: " << rms << std::endl;
  oops::Log::debug() << "rms from config: " << rms_conf << std::endl;

  BOOST_CHECK_CLOSE(rms, rms_conf, tol);
}

// -----------------------------------------------------------------------------

template <typename MODEL> class GeometryIterator : public oops::Test {
 public:
  GeometryIterator() {}
  virtual ~GeometryIterator() {}
 private:
  std::string testid() const {return "test::GeometryIterator<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/GeometryIterator");

    ts->add(BOOST_TEST_CASE(&testConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testIterator<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_GEOMETRYITERATOR_H_
