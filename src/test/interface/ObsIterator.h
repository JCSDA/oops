/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_INTERFACE_OBSITERATOR_H_
#define TEST_INTERFACE_OBSITERATOR_H_

#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/Configuration.h"
#include "eckit/testing/Test.h"

#include "oops/interface/GeometryIterator.h"
#include "oops/interface/ObsSpace.h"
#include "oops/runs/Test.h"

#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------
/*! \brief Tests of ObsSpace::begin/end; ObsIterator ctor, ==/!=/++ operators
 *         and dereferencing operator.
 *
 * \details testBasic tests the following:
 *
 * 1. Initialize ObsIterator to ObsSpace::begin() and check equality
 * 2. Initialize ObsIterator to ObsSpace::end() and check equality
 * 3. Check inequality of the two iterators
 * 4. Print out the begin iterator, to "test" print method
 * 5. Loop over iterator (test operator++), compute number of iterations, compare
 *    with reference
 * 6. Compare the first two Points with reference.
 */
template <typename OBS> void testBasic() {
  typedef oops::GeometryIterator<OBS>   ObsIterator_;
  typedef ObsTestsFixture<OBS>          Test_;

  for (size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    const double tol = Test_::config(jj).getDouble("obs iterator test.tolerance");

    // Initialize two iterators (begin and end), test equality
    ObsIterator_ iter1 = Test_::obspace()[jj].begin();
    EXPECT(iter1 == Test_::obspace()[jj].begin());

    ObsIterator_ iter2 = Test_::obspace()[jj].end();
    EXPECT(iter2 == Test_::obspace()[jj].end());

    // For all current use cases begin() shouldn't be the same as end(); test it
    EXPECT(iter1 != iter2);

    // At least test that nothing fails on print
    oops::Log::info() << "ObsSpace::begin " << iter1 << std::endl;

    // loop over all points in iterator, compute number of those points, and compare
    // with reference
    size_t nlocs = 0;
    for (ObsIterator_ ii = Test_::obspace()[jj].begin(); ii != Test_::obspace()[jj].end(); ++ii) {
      nlocs++;
    }
    size_t nlocs_ref = Test_::config(jj).getUnsigned("obs iterator test.reference nlocs");
    EXPECT_EQUAL(nlocs, nlocs_ref);

    // test that the begin() point is the same as reference
    double lon1 = Test_::config(jj).getDouble("obs iterator test.lon1");
    double lat1 = Test_::config(jj).getDouble("obs iterator test.lat1");
    const eckit::geometry::Point3 point1(lon1, lat1, 0.0);

    EXPECT((*iter1).distance(point1) <= tol);

    // test that the point after begin() is the same as reference
    ++iter1;
    double lon2 = Test_::config(jj).getDouble("obs iterator test.lon2");
    double lat2 = Test_::config(jj).getDouble("obs iterator test.lat2");
    const eckit::geometry::Point3 point2(lon2, lat2, 0.0);
    EXPECT((*iter1).distance(point2) <= tol);
  }
}


// -----------------------------------------------------------------------------

template <typename OBS> class ObsIterator : public oops::Test {
  typedef ObsTestsFixture<OBS> Test_;
 public:
  using oops::Test::Test;
  virtual ~ObsIterator() = default;

 private:
  std::string testid() const override {return "test::ObsIterator<" + OBS::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/ObsIterator/testBasic")
      { testBasic<OBS>(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSITERATOR_H_
