/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_MODELAUXCONTROL_H_
#define TEST_INTERFACE_MODELAUXCONTROL_H_

#include <string>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/runs/Test.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/ModelAuxControl.h"
#include "test/TestEnvironment.h"
#include "eckit/config/LocalConfiguration.h"

namespace test {

// -----------------------------------------------------------------------------
template <typename MODEL> class ModelAuxControlFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL>        Geometry_;
  typedef oops::ModelAuxControl<MODEL> ModelAux_;

 public:
  static const eckit::Configuration & config() {return *getInstance().conf_;}
  static const Geometry_    & resol()  {return *getInstance().resol_;}

 private:
  static ModelAuxControlFixture<MODEL>& getInstance() {
    static ModelAuxControlFixture<MODEL> theModelAuxControlFixture;
    return theModelAuxControlFixture;
  }

  ModelAuxControlFixture() {
    conf_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "ModelBias"));

    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "Geometry");
    resol_.reset(new Geometry_(resolConfig));
  }

  ~ModelAuxControlFixture() {}

  boost::scoped_ptr<const eckit::LocalConfiguration>  conf_;
  boost::scoped_ptr<Geometry_>     resol_;
};

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef ModelAuxControlFixture<MODEL>   Test_;
  typedef oops::ModelAuxControl<MODEL>    ModelAux_;

  boost::scoped_ptr<ModelAux_> bias(new ModelAux_(Test_::resol(), Test_::config()));
  BOOST_CHECK(bias.get());

  bias.reset();
  BOOST_CHECK(!bias.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testCopyConstructor() {
  typedef ModelAuxControlFixture<MODEL>   Test_;
  typedef oops::ModelAuxControl<MODEL>    ModelAux_;

  boost::scoped_ptr<ModelAux_> bias(new ModelAux_(Test_::resol(), Test_::config()));

  boost::scoped_ptr<ModelAux_> other(new ModelAux_(*bias));
  BOOST_CHECK(other.get());

  other.reset();
  BOOST_CHECK(!other.get());

  BOOST_CHECK(bias.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testChangeRes() {
  typedef ModelAuxControlFixture<MODEL>   Test_;
  typedef oops::ModelAuxControl<MODEL>    ModelAux_;

  boost::scoped_ptr<ModelAux_> bias(new ModelAux_(Test_::resol(), Test_::config()));

  boost::scoped_ptr<ModelAux_> other(new ModelAux_(Test_::resol(), *bias));
  BOOST_CHECK(other.get());

  other.reset();
  BOOST_CHECK(!other.get());

  BOOST_CHECK(bias.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> class ModelAuxControl : public oops::Test {
 public:
  ModelAuxControl() {}
  virtual ~ModelAuxControl() {}
 private:
  std::string testid() const {return "test::ModelAuxControl<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/ModelAuxControl");

    ts->add(BOOST_TEST_CASE(&testConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testCopyConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testChangeRes<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_MODELAUXCONTROL_H_
