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

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/shared_ptr.hpp>

#include "eckit/testing/Test.h"
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
      eckit::LocalConfiguration osconf(conf[jj], "ObsSpace");
      boost::shared_ptr<ObsSpace_> tmp(new ObsSpace_(osconf, bgn, end));
      ospaces_.push_back(tmp);
      eckit::LocalConfiguration ObsDataInConf;
      osconf.get("ObsDataIn", ObsDataInConf);
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
    std::unique_ptr<ObsVector_> ov(new ObsVector_(*Test_::obspace()[jj]));
    EXPECT(ov.get());

    ov.reset();
    EXPECT(!ov.get());
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testCopyConstructor() {
  typedef ObsVecFixture<MODEL>  Test_;
  typedef oops::ObsVector<MODEL>  ObsVector_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    std::unique_ptr<ObsVector_> ov(new ObsVector_(*Test_::obspace()[jj]));
    std::unique_ptr<ObsVector_> other(new ObsVector_(*ov));
    EXPECT(other.get());

    other.reset();
    EXPECT(!other.get());

    EXPECT(ov.get());
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
    EXPECT(ovov2 > zero);

    ov.zero();

    const double zz = dot_product(ov, ov);
    EXPECT(zz == zero);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsVector : public oops::Test {
 public:
  ObsVector() {}
  virtual ~ObsVector() {}
 private:
  std::string testid() const {return "test::ObsVector<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/ObsVector/testConstructor")
      { testConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/ObsVector/testCopyConstructor")
      { testCopyConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/ObsVector/testNotZero")
      { testNotZero<MODEL>(); });
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSVECTOR_H_
