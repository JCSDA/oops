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

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/Configuration.h"
#include "eckit/testing/Test.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/LocalIncrement.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"

#include "test/interface/GeometryFixture.h"
#include "test/interface/Increment.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------
/*! \brief Tests of Geometry::begin/end; GeometryIterator ctor and ==/!= operators
 *
 * \details testBasic tests the following:
 *
 * 1. Initialize GeometryIterator to Geometry::begin() and check equality
 * 2. Initialize GeometryIterator to Geometry::end() and check equality
 * 3. Check inequality of the two iterators
 * 4. Print out the begin iterator, to "test" print method
 */
template <typename MODEL> void testBasic() {
  typedef oops::GeometryIterator<MODEL>   GeometryIterator_;
  typedef oops::Geometry<MODEL>           Geometry_;

  Geometry_ geom(GeometryFixture<MODEL>::getParameters(), oops::mpi::world(), oops::mpi::myself());

  GeometryIterator_ iter1 = geom.begin();
  EXPECT(iter1 == geom.begin());

  GeometryIterator_ iter2 = geom.end();
  EXPECT(iter2 == geom.end());

  // For all current use cases begin() shouldn't be the same as end(); test it
  EXPECT(iter1 != iter2);

  // At least test that nothing fails on print
  oops::Log::test() << "Geometry::begin " << iter1 << std::endl;
}


// -----------------------------------------------------------------------------
/*! \brief Test of GeometryIterator::operator++, Increment::getLocal and Increment::setLocal
 *
 * \details testGetSetLocal tests the following:
 *
 * Initialize dx1 to non-zero random values, initialize dx2 to zero.
 * Loop through gridpoints using GeometryIterator operator++, assign dx2 to dx1 gridpoint
 * by gridpoint using get/setLocal. Check that the two increments are the same
 * in the end.
 */
template <typename MODEL> void testGetSetLocal() {
  typedef oops::Geometry<MODEL>          Geometry_;
  typedef oops::GeometryIterator<MODEL>  GeometryIterator_;
  typedef oops::Increment<MODEL>         Increment_;
  typedef IncrementFixture<MODEL>        Test_;

  Geometry_ geom(GeometryFixture<MODEL>::getParameters(), oops::mpi::world(), oops::mpi::myself());

  // randomize increment dx1
  Increment_ dx1(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx1.random();
  EXPECT(dx1.norm() != 0.0);
  oops::Log::test() << "Increment dx1 (random): " << dx1 << std::endl;
  // zero out increment dx2
  Increment_ dx2(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx2.zero();
  EXPECT(dx2.norm() == 0.0);
  oops::Log::test() << "Increment dx2 (zero): " << dx2 << std::endl;

  for (GeometryIterator_ i = geom.begin(); i != geom.end(); ++i) {
    // get value for i-th gridpoint from dx1
    oops::LocalIncrement gp = dx1.getLocal(i);
    oops::Log::debug() << *i << gp << std::endl;
    // set this value for i-th gridpoint in dx2
    dx2.setLocal(gp, i);
  }
  oops::Log::test() << "Increment dx2 after dx2=dx1 (at every point): " << dx2 << std::endl;
  EXPECT(dx2.norm() != 0.0);
  // compare two increments
  dx2 -= dx1;
  EXPECT(dx2.norm() == 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> class GeometryIterator : public oops::Test {
 public:
  GeometryIterator() {}
  virtual ~GeometryIterator() {IncrementFixture<MODEL>::reset();}

 private:
  std::string testid() const override {return "test::GeometryIterator<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/GeometryIterator/testBasic")
      { testBasic<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/testGetSetLocal")
      { testGetSetLocal<MODEL>(); });
  }

  void clear() const override {}
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_GEOMETRYITERATOR_H_
