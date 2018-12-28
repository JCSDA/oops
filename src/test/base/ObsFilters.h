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

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/ObsFilterBase.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/Locations.h"
#include "oops/interface/ObsAuxControl.h"
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
  typedef oops::Locations<MODEL>         Locations_;
  typedef oops::ObsAuxControl<MODEL>     ObsAuxCtrl_;
  typedef oops::ObsFilterBase<MODEL>     ObsFilterBase_;
  typedef oops::ObsOperator<MODEL>       ObsOperator_;
  typedef oops::ObsVector<MODEL>         ObsVector_;

  const eckit::LocalConfiguration obsconf(TestEnvironment::config(), "Observations");
  std::vector<eckit::LocalConfiguration> typeconfs;
  obsconf.get("ObsTypes", typeconfs);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
//  Would not need hop for this test if using precomputed geovals
    ObsOperator_ hop(Test_::obspace()[jj], typeconfs[jj]);
    eckit::LocalConfiguration gconf(typeconfs[jj], "GeoVaLs");
    Locations_ locs(hop.locations(Test_::tbgn(), Test_::tend()));
    const GeoVaLs_ gval(gconf, hop.variables());

    eckit::LocalConfiguration biasConf;
    typeconfs[jj].get("ObsBias", biasConf);
    const ObsAuxCtrl_ ybias(biasConf);

    ObsVector_ ovec(Test_::obspace()[jj], hop.observed());

//  Prepare storage for QC flags using PreQC filter
    const std::string qcname = "TestQC";
    eckit::LocalConfiguration preconf;
    preconf.set("Filter", "PreQC");
    preconf.set("observed", hop.observed().variables());
    preconf.set("QCname", qcname);
    boost::shared_ptr<ObsFilterBase_> preqc(  // All work is done in constructor
      oops::FilterFactory<MODEL>::create(Test_::obspace()[jj], preconf));

//  Get filters configurations
    std::vector<eckit::LocalConfiguration> filtconf;
    typeconfs[jj].get("ObsFilters", filtconf);
    oops::Log::debug() << "test filt conf " << filtconf[jj] << std::endl;

//  Test filters
    for (std::size_t jf = 0; jf < filtconf.size(); ++jf) {
      filtconf[jf].set("QCname", qcname);
      boost::shared_ptr<ObsFilterBase_> filter(
        oops::FilterFactory<MODEL>::create(Test_::obspace()[jj], filtconf[jf]));

      filter->priorFilter(gval);
      hop.simulateObs(gval, ovec, ybias);
      filter->postFilter(ovec);
    }

    const double tol = typeconfs[jj].getDouble("tolerance");
    if (typeconfs[jj].has("vecequiv")) {
      // if reference h(x) is saved in file as a vector, read from file
      // and compare the norm of difference to zero
      ObsVector_ ovec_ref(ovec, false);
      ovec_ref.read(typeconfs[jj].getString("vecequiv"));
      ovec_ref -= ovec;
      const double zz = ovec_ref.rms();
      oops::Log::info() << "Vector difference between reference and computed: " <<
                           ovec_ref;
      BOOST_CHECK_SMALL(zz, tol);
    } else {
      // else compare h(x) norm to the norm from the config
      const double zz = ovec.rms();
      const double xx = typeconfs[jj].getDouble("rmsequiv");
      BOOST_CHECK_CLOSE(xx, zz, tol);
    }

  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> class ObsFilters : public oops::Test {
 public:
  ObsFilters() {}
  virtual ~ObsFilters() {}
 private:
  std::string testid() const {return "test::ObsFilters<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("base/ObsFilters");

    ts->add(BOOST_TEST_CASE(&testFilters<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_BASE_OBSFILTERS_H_
