/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_MODELAUXCOVARIANCE_H_
#define TEST_INTERFACE_MODELAUXCOVARIANCE_H_

#include <string>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/runs/Test.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/ModelAuxCovariance.h"
#include "test/TestEnvironment.h"
#include "eckit/config/LocalConfiguration.h"

namespace test {

// -----------------------------------------------------------------------------
template <typename MODEL> class ModelAuxCovarianceFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL>           Geometry_;
  typedef oops::ModelAuxCovariance<MODEL> Covariance_;

 public:
  static const eckit::Configuration & config() {return *getInstance().conf_;}
  static const Geometry_    & resol()  {return *getInstance().resol_;}

 private:
  static ModelAuxCovarianceFixture<MODEL>& getInstance() {
    static ModelAuxCovarianceFixture<MODEL> theModelAuxCovarianceFixture;
    return theModelAuxCovarianceFixture;
  }

  ModelAuxCovarianceFixture() {
    conf_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "ModelBiasCovariance"));

    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "Geometry");
    resol_.reset(new Geometry_(resolConfig));
  }

  ~ModelAuxCovarianceFixture() {}

  boost::scoped_ptr<const eckit::LocalConfiguration>  conf_;
  boost::scoped_ptr<Geometry_>     resol_;
};

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef ModelAuxCovarianceFixture<MODEL>   Test_;
  typedef oops::ModelAuxCovariance<MODEL>    Covariance_;

  boost::scoped_ptr<Covariance_> bias(new Covariance_(Test_::config(), Test_::resol()));
  BOOST_CHECK(bias.get());

  bias.reset();
  BOOST_CHECK(!bias.get());
}

// -----------------------------------------------------------------------------
//  void linearize(const ModelAuxControl_ &, const Geometry_ &);
//  void multiply(const ModelAuxIncrement_ &, ModelAuxIncrement_ &) const;
//  void invMult(const ModelAuxIncrement_ &, ModelAuxIncrement_ &) const;
//  void randomize(ModelAuxIncrement_ &) const;
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

template <typename MODEL> class ModelAuxCovariance : public oops::Test {
 public:
  ModelAuxCovariance() {}
  virtual ~ModelAuxCovariance() {}
 private:
  std::string testid() const {return "test::ModelAuxCovariance<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/ModelAuxCovariance");

    ts->add(BOOST_TEST_CASE(&testConstructor<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_MODELAUXCOVARIANCE_H_
