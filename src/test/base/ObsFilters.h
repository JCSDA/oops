/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_BASE_OBSFILTERS_H_
#define TEST_BASE_OBSFILTERS_H_

#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/ObsFilters.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsVector.h"
#include "oops/runs/Test.h"
#include "oops/util/Logger.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL> void testFilters() {
  typedef ObsTestsFixture<MODEL> Test_;
  typedef oops::GeoVaLs<MODEL>           GeoVaLs_;
  typedef oops::ObsAuxControl<MODEL>     ObsAuxCtrl_;
  typedef oops::ObsFilters<MODEL>        ObsFilters_;
  typedef oops::ObsOperator<MODEL>       ObsOperator_;
  typedef oops::ObsVector<MODEL>         ObsVector_;

  const eckit::LocalConfiguration obsconf(TestEnvironment::config(), "Observations");
  std::vector<eckit::LocalConfiguration> typeconfs;
  obsconf.get("ObsTypes", typeconfs);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    eckit::LocalConfiguration obsopconf(typeconfs[jj], "ObsOperator");
    ObsOperator_ hop(Test_::obspace()[jj], obsopconf);

    eckit::LocalConfiguration gconf(typeconfs[jj], "GeoVaLs");
    const GeoVaLs_ gval(gconf, hop.variables());

    const ObsAuxCtrl_ ybias(typeconfs[jj]);

//  Allocate memory for tests
    ObsVector_ ovec(Test_::obspace()[jj]);
    boost::shared_ptr<oops::ObsDataVector<MODEL, float> > obserr
      (new oops::ObsDataVector<MODEL, float>(Test_::obspace()[jj],
               Test_::obspace()[jj].obsvariables(), "ObsError"));
    boost::shared_ptr<oops::ObsDataVector<MODEL, int> >
      qcflags(new oops::ObsDataVector<MODEL, int>  (Test_::obspace()[jj],
               Test_::obspace()[jj].obsvariables()));

//  Create filters
    ObsFilters_ filters(Test_::obspace()[jj], typeconfs[jj], qcflags, obserr);

//  Run filters
    filters.preProcess();
    filters.priorFilter(gval);
    hop.simulateObs(gval, ovec, ybias);
    filters.postFilter(ovec);

//  Compare with known results
    if (typeconfs[jj].has("qcBenchmark")) {
      const std::string qcBenchmarkName = typeconfs[jj].getString("qcBenchmark");

      oops::ObsDataVector<MODEL, int> qcBenchmark(Test_::obspace()[jj],
                           Test_::obspace()[jj].obsvariables(), qcBenchmarkName);

      bool same = compareFlags(*qcflags, qcBenchmark);
      EXPECT(same);
    } else {
      const int passedBenchmark = typeconfs[jj].getInt("passedBenchmark");
      const int passed = numZero(*qcflags);
      EXPECT(passed == passedBenchmark);
    }
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsFilters : public oops::Test {
 public:
  ObsFilters() {}
  virtual ~ObsFilters() {}
 private:
  std::string testid() const {return "test::ObsFilters<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("base/ObsFilters/testFilters")
      { testFilters<MODEL>(); });
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_BASE_OBSFILTERS_H_
