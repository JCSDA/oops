/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_LINEARMODEL_H_
#define TEST_INTERFACE_LINEARMODEL_H_

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include <cmath>
#include <vector>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "oops/runs/Test.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/TrajectorySaver.h"
#include "oops/generic/instantiateTlmFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/LinearModel.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/interface/State.h"
#include "oops/interface/Variables.h"
#include "test/TestEnvironment.h"
#include "eckit/config/LocalConfiguration.h"
#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/dot_product.h"


using oops::Log;

namespace test {

// =============================================================================

template <typename MODEL> class LinearModelFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL>          Geometry_;
  typedef oops::Increment<MODEL>         Increment_;
  typedef oops::LinearModel<MODEL>       LinearModel_;
  typedef oops::Model<MODEL>             Model_;
  typedef oops::ModelAuxControl<MODEL>   ModelAux_;
  typedef oops::ModelAuxIncrement<MODEL> ModelAuxIncr_;
  typedef oops::State<MODEL>             State_;
  typedef oops::Variables<MODEL>         Variables_;
  typedef oops::ModelSpaceCovarianceBase<MODEL> Covariance_;

 public:
  static const eckit::Configuration & test()   {return *getInstance().test_;}
  static const Geometry_        & resol()      {return *getInstance().resol_;}
  static const Variables_       & ctlvars()    {return *getInstance().ctlvars_;}
  static const util::DateTime   & time()       {return *getInstance().time_;}
  static const Covariance_      & covariance() {return *getInstance().B_;}
  static const Model_           & model()      {return *getInstance().model_;}
  static const State_           & xref()       {return *getInstance().xref_;}
  static const ModelAux_        & bias()       {return *getInstance().bias_;}
  static const ModelAuxIncr_    & dbias()      {return *getInstance().dbias_;}
  static const LinearModel_     & tlm()        {return *getInstance().tlm_;}

 private:
  static LinearModelFixture<MODEL>& getInstance() {
    static LinearModelFixture<MODEL> theLinearModelFixture;
    return theLinearModelFixture;
  }

  LinearModelFixture<MODEL>() {
    test_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "LinearModelTest"));
    const util::Duration len(test_->getString("fclength"));

    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "Geometry");
    resol_.reset(new Geometry_(resolConfig));

    const eckit::LocalConfiguration varConfig(TestEnvironment::config(), "Variables");
    ctlvars_.reset(new Variables_(varConfig));

    const eckit::LocalConfiguration iniConf(TestEnvironment::config(), "State");
    xref_.reset(new State_(*resol_, iniConf));
    time_.reset(new util::DateTime(xref_->validTime()));

    const eckit::LocalConfiguration biasConf(TestEnvironment::config(), "ModelBias");
    bias_.reset(new ModelAux_(*resol_, biasConf));
    dbias_.reset(new ModelAuxIncr_(*resol_, biasConf));

    const eckit::LocalConfiguration nlConf(TestEnvironment::config(), "Model");
    model_.reset(new Model_(*resol_, nlConf));

//  Create a covariance matrix
    oops::instantiateCovarFactory<MODEL>();
    const eckit::LocalConfiguration covar(TestEnvironment::config(), "Covariance");
    Covariance_ * Bptr = oops::CovarianceFactory<MODEL>::create(covar, *resol_, *ctlvars_, *xref_);
    Bptr->linearize(*xref_, *resol_);
    B_.reset(Bptr);

//  Linear model configuration
    tlConf_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "LinearModel"));

//  Setup trajectory for TL and AD
    oops::instantiateTlmFactory<MODEL>();
    boost::ptr_vector<LinearModel_> tlmvec;
    oops::PostProcessor<State_> post;
    post.enrollProcessor(new oops::TrajectorySaver<MODEL>(*xref_, *tlConf_, *resol_, *bias_, tlmvec));
    State_ xx(*xref_);
    model_->forecast(xx, *bias_, len, post);
    tlm_.reset(tlmvec.release(tlmvec.begin()).release());
  }

  ~LinearModelFixture<MODEL>() {}

  boost::scoped_ptr<const eckit::LocalConfiguration>   test_;
  boost::scoped_ptr<const eckit::LocalConfiguration>   tlConf_;
  boost::scoped_ptr<const Geometry_>      resol_;
  boost::scoped_ptr<const util::DateTime> time_;
  boost::scoped_ptr<const Variables_>     ctlvars_;
  boost::scoped_ptr<const State_>         xref_;
  boost::scoped_ptr<const Model_>         model_;
  boost::scoped_ptr<const ModelAux_>      bias_;
  boost::scoped_ptr<const ModelAuxIncr_>  dbias_;
  boost::scoped_ptr<const Covariance_>    B_;
  boost::scoped_ptr<const LinearModel_>   tlm_;
};

// =============================================================================

template <typename MODEL> void testLinearModelConstructor() {
  typedef LinearModelFixture<MODEL>   Test_;

  const util::Duration zero(0);
  BOOST_CHECK(Test_::tlm().timeResolution() > zero);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearModelZeroLength() {
  typedef LinearModelFixture<MODEL>      Test_;
  typedef oops::Increment<MODEL>         Increment_;
  typedef oops::ModelAuxIncrement<MODEL> ModelAuxIncr_;

  const util::DateTime vt(Test_::time());
  const util::Duration zero(0);

  Increment_ dx(Test_::resol(), Test_::ctlvars(), vt);
  Test_::covariance().randomize(dx);
  ModelAuxIncr_ daux(Test_::dbias());
  const double ininorm = dx.norm();
  BOOST_CHECK(ininorm > 0.0);

  Test_::tlm().forecastTL(dx, daux, zero);
  BOOST_CHECK_EQUAL(dx.validTime(), vt);
  BOOST_CHECK_EQUAL(dx.norm(), ininorm);

  Test_::tlm().forecastAD(dx, daux, zero);
  BOOST_CHECK_EQUAL(dx.validTime(), vt);
  BOOST_CHECK_EQUAL(dx.norm(), ininorm);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearModelZeroPert() {
  typedef LinearModelFixture<MODEL>      Test_;
  typedef oops::Increment<MODEL>         Increment_;
  typedef oops::ModelAuxIncrement<MODEL> ModelAuxIncr_;

  const util::Duration len(Test_::test().getString("fclength"));
  const util::DateTime t1(Test_::time());
  const util::DateTime t2(t1 + len);
  BOOST_CHECK(t2 > t1);

  Increment_ dx(Test_::resol(), Test_::ctlvars(), t1);
  ModelAuxIncr_ daux(Test_::dbias());

  dx.zero();
  daux.zero();
  BOOST_CHECK_EQUAL(dx.norm(), 0.0);
  Test_::tlm().forecastTL(dx, daux, len);
  BOOST_CHECK_EQUAL(dx.validTime(), t2);
  BOOST_CHECK_EQUAL(dx.norm(), 0.0);

  dx.zero();
  daux.zero();
  BOOST_CHECK_EQUAL(dx.norm(), 0.0);
  Test_::tlm().forecastAD(dx, daux, len);
  BOOST_CHECK_EQUAL(dx.validTime(), t1);
  BOOST_CHECK_EQUAL(dx.norm(), 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearModelLinearity() {
  typedef LinearModelFixture<MODEL>      Test_;
  typedef oops::Increment<MODEL>         Increment_;
  typedef oops::ModelAuxIncrement<MODEL> ModelAuxIncr_;

  const util::Duration len(Test_::test().getString("fclength"));
  const util::DateTime t1(Test_::time());
  const util::DateTime t2(t1 + len);
  BOOST_CHECK(t2 > t1);
  const double zz = 3.1415;

  Increment_ dx1(Test_::resol(), Test_::ctlvars(), t1);
  Test_::covariance().randomize(dx1);
  ModelAuxIncr_ daux1(Test_::dbias());
  BOOST_CHECK(dx1.norm() > 0.0);

  Increment_ dx2(dx1);
  ModelAuxIncr_ daux2(daux1);

  Test_::tlm().forecastTL(dx1, daux1, len);
  BOOST_CHECK_EQUAL(dx1.validTime(), t2);
  dx1 *= zz;
  daux1 *= zz;

  dx2 *= zz;
  daux2 *= zz;
  Test_::tlm().forecastTL(dx2, daux2, len);
  BOOST_CHECK_EQUAL(dx2.validTime(), t2);

  const double tol = Test_::test().getDouble("toleranceAD");
  BOOST_CHECK_CLOSE(dx1.norm(), dx2.norm(), tol);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearApproximation() {
  typedef LinearModelFixture<MODEL>      Test_;
  typedef oops::Increment<MODEL>         Increment_;
  typedef oops::State<MODEL>             State_;

  const util::Duration len(Test_::test().getString("fclength"));
  const util::DateTime t1(Test_::time());
  const util::DateTime t2(t1 + len);
  BOOST_CHECK(t2 > t1);

  Increment_ dx0(Test_::resol(), Test_::ctlvars(), t1);
  Test_::covariance().randomize(dx0);
  BOOST_CHECK(dx0.norm() > 0.0);

  Increment_ dx(dx0);
  Test_::tlm().forecastTL(dx, Test_::dbias(), len);
  const double dxnorm = dx.norm();

  oops::PostProcessor<State_> post;
  State_ xx0(Test_::xref());
  Test_::model().forecast(xx0, Test_::bias(), len, post);

  const unsigned int ntest = Test_::test().getInt("testiterTL");
  double zz = 1.0;
  std::vector<double> errors;
  for (unsigned int jtest = 0; jtest < ntest; ++jtest) {
    State_ xx(Test_::xref());
    Increment_ pert(dx0);
    pert *= zz;
    xx += pert;
    Test_::model().forecast(xx, Test_::bias(), len, post);

    Increment_ diff(Test_::resol(), Test_::ctlvars(), t2);
    diff.diff(xx, xx0);
    const double difnorm = diff.norm();
    const double err = zz * dxnorm / difnorm;
    Increment_ derr(dx);
    derr *= zz;
    derr -= diff;
    const double errnorm = derr.norm();
    errors.push_back(errnorm / difnorm);
    Log::test() << "TL error = " << std::setprecision(16) << err
                << ", relative error = " << errnorm / difnorm << std::endl;
    zz /= 10.0;
  }

// Analyze results
  const double approx = *std::min_element(errors.begin(), errors.end());
  Log::test() << "Test TL min error = " << approx << std::endl;
  const double tol = Test_::test().getDouble("toleranceTL");
  BOOST_CHECK(approx < tol);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearModelAdjoint() {
  typedef LinearModelFixture<MODEL>      Test_;
  typedef oops::Increment<MODEL>         Increment_;
  typedef oops::ModelAuxIncrement<MODEL> ModelAuxIncr_;

  const util::Duration len(Test_::test().getString("fclength"));
  const util::DateTime t1(Test_::time());
  const util::DateTime t2(t1 + len);
  BOOST_CHECK(t2 > t1);

  Increment_ dx11(Test_::resol(), Test_::ctlvars(), t1);
  Test_::covariance().randomize(dx11);
  ModelAuxIncr_ daux1(Test_::dbias());
  BOOST_CHECK(dx11.norm() > 0.0);
  Increment_ dx12(dx11);
  Test_::tlm().forecastTL(dx12, daux1, len);
  BOOST_CHECK(dx12.norm() > 0.0);

  Increment_ dx22(Test_::resol(), Test_::ctlvars(), t2);
  Test_::covariance().randomize(dx22);
  ModelAuxIncr_ daux2(Test_::dbias());
  BOOST_CHECK(dx22.norm() > 0.0);
  Increment_ dx21(dx22);
  Test_::tlm().forecastAD(dx21, daux2, len);
  BOOST_CHECK(dx21.norm() > 0.0);

  BOOST_CHECK(dx11.norm() != dx22.norm());
  BOOST_CHECK_EQUAL(dx11.validTime(), t1);
  BOOST_CHECK_EQUAL(dx21.validTime(), t1);
  BOOST_CHECK_EQUAL(dx12.validTime(), t2);
  BOOST_CHECK_EQUAL(dx22.validTime(), t2);

  const double dot1 = dot_product(dx11, dx21);
  const double dot2 = dot_product(dx12, dx22);
  const double tol = Test_::test().getDouble("toleranceAD");
  BOOST_CHECK_CLOSE(dot1, dot2, tol);
}

// =============================================================================

template <typename MODEL> class LinearModel : public oops::Test {
 public:
  LinearModel() {}
  virtual ~LinearModel() {}
 private:
  std::string testid() const {return "test::LinearModel<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/LinearModel");

    ts->add(BOOST_TEST_CASE(&testLinearModelConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testLinearModelZeroLength<MODEL>));
    ts->add(BOOST_TEST_CASE(&testLinearModelZeroPert<MODEL>));
    ts->add(BOOST_TEST_CASE(&testLinearModelLinearity<MODEL>));
    ts->add(BOOST_TEST_CASE(&testLinearApproximation<MODEL>));
    ts->add(BOOST_TEST_CASE(&testLinearModelAdjoint<MODEL>));
    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_LINEARMODEL_H_
