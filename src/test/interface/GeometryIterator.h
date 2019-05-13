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

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "eckit/testing/Test.h"
#include "oops/base/GridPoint.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/interface/State.h"
#include "oops/runs/Test.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef oops::GeometryIterator<MODEL>   GeometryIterator_;
  typedef oops::Geometry<MODEL>           Geometry_;

  const eckit::LocalConfiguration
       geomConfig(TestEnvironment::config(), "Geometry");
  Geometry_ geom(geomConfig);

  boost::scoped_ptr<GeometryIterator_> iter(new GeometryIterator_(geom.begin()));
  EXPECT(iter.get());

  iter.reset();
  EXPECT(!iter.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testIterator() {
  typedef oops::Geometry<MODEL>          Geometry_;
  typedef oops::GeometryIterator<MODEL>  GeometryIterator_;
  typedef oops::State<MODEL>             State_;

  const eckit::LocalConfiguration
       geomConfig(TestEnvironment::config(), "Geometry");
  const eckit::LocalConfiguration
       stateConfig(TestEnvironment::config(), "State");

  const oops::Variables vars(stateConfig);

  Geometry_ geom(geomConfig);
  State_ state(geom, vars, stateConfig);

  const eckit::LocalConfiguration
        iterConfig(TestEnvironment::config(), "GeometryIterator");
  const double rms_conf = iterConfig.getDouble("rms");
  const double tol = iterConfig.getDouble("tolerance");

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

  oops::Log::debug() << n << " rms from iterator: " << std::fixed << std::setprecision(8)
                     << rms << std::endl;
  oops::Log::debug() << "rms from config: " << std::fixed << std::setprecision(8)
                     << rms_conf << std::endl;

  EXPECT(oops::is_close(rms, rms_conf, tol));
}

// -----------------------------------------------------------------------------

template <typename MODEL> class GeometryIterator : public oops::Test {
 public:
  GeometryIterator() {}
  virtual ~GeometryIterator() {}

 private:
  std::string testid() const {return "test::GeometryIterator<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/GeometryIterator/testConstructor")
      { testConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/GeometryIterator/testIterator")
      { testIterator<MODEL>(); });
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_GEOMETRYITERATOR_H_
