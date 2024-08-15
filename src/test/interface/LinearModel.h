/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Crown Copyright 2024, the Met Office.
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
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/Increment4D.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/LinearModel.h"
#include "oops/base/Model.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/State4D.h"
#include "oops/base/TrajectorySaver.h"
#include "oops/base/Variables.h"
#include "oops/generic/instantiateLinearModelFactory.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Duration.h"
#include "test/TestEnvironment.h"

namespace test {

// =============================================================================

template <typename MODEL> class LinearModelFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL>          Geometry_;
  typedef oops::Increment<MODEL>         Increment_;
  typedef oops::Increment4D<MODEL>       Increment4D_;
  typedef oops::LinearModel<MODEL>       LinearModel_;
  typedef oops::Model<MODEL>             Model_;
  typedef oops::ModelAuxControl<MODEL>   ModelAux_;
  typedef oops::ModelAuxIncrement<MODEL> ModelAuxIncr_;
  typedef oops::State<MODEL>             State_;
  typedef oops::State4D<MODEL>           State4D_;
  typedef oops::ModelSpaceCovarianceBase<MODEL> Covariance_;

 public:
  static const eckit::Configuration & test()   {return *getInstance().test_;}
  static const Geometry_        & resol()      {return *getInstance().resol_;}
  static const oops::Variables  & ctlvars()    {return *getInstance().ctlvars_;}
  static const util::DateTime   & time()       {return *getInstance().time_;}
  static const Covariance_      & covariance() {return *getInstance().B_;}
  static const Model_           & model()      {return *getInstance().model_;}
  static const State4D_         & xref()       {return *getInstance().xref_;}
  static const ModelAux_        & bias()       {return *getInstance().bias_;}
  static const ModelAuxIncr_    & dbias()      {return *getInstance().dbias_;}
  static const LinearModel_     & tlm()        {return *getInstance().tlm_;}
  static void reset() {
    getInstance().tlm_.reset();
    getInstance().B_.reset();
    getInstance().time_.reset();
    getInstance().xref_.reset();
    getInstance().model_.reset();
    getInstance().dbias_.reset();
    getInstance().bias_.reset();
    getInstance().ctlvars_.reset();
    getInstance().resol_.reset();
    getInstance().test_.reset();
  }

  /// \brief Returns the instance of the LinearModelFixture.
  /// \detail The first call requires a pointer to the MPI communicator,
  ///         in order to initialise the Geometry.
  static LinearModelFixture<MODEL>& getInstance(
      const eckit::mpi::Comm * const comm_setter = nullptr) {
    // Create instance using the communicator.
    static LinearModelFixture<MODEL> theLinearModelFixture([=]() -> const eckit::mpi::Comm & {
      ASSERT_MSG(comm_setter != nullptr,
                 "LinearModelFixture::getInstance: ERROR: MPI communicator required on"
                 " the first call to `getInstance`.");
      return *comm_setter;
    }());

    return theLinearModelFixture;
  }

 private:
  explicit LinearModelFixture<MODEL>(const eckit::mpi::Comm & comm) {
    test_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "linear model test"));
    const util::Duration len(test_->getString("forecast length"));

    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "geometry");
    resol_.reset(new Geometry_(resolConfig, comm));

    ctlvars_.reset(new oops::Variables(TestEnvironment::config(), "analysis variables"));

    const eckit::LocalConfiguration biasConf(TestEnvironment::config(), "model aux control");
    bias_.reset(new ModelAux_(*resol_, biasConf));

    const eckit::LocalConfiguration dbiasConf =
          TestEnvironment::config().getSubConfiguration("model aux error");
    dbias_.reset(new ModelAuxIncr_(*resol_, dbiasConf));

    const eckit::LocalConfiguration nlConf(TestEnvironment::config(), "model");
    model_.reset(new Model_(*resol_, nlConf));

    const eckit::LocalConfiguration iniConf(TestEnvironment::config(), "initial condition");
    xref_.reset(new State4D_(*resol_, iniConf));
    time_.reset(new util::DateTime(xref_->times()[0]));

//  Create a covariance matrix
    oops::instantiateCovarFactory<MODEL>();
    const eckit::LocalConfiguration covar(TestEnvironment::config(), "background error");
    B_.reset(oops::CovarianceFactory<MODEL>::create(*resol_, *ctlvars_, covar, *xref_, *xref_));

//  Linear model configuration
    tlConf_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "linear model"));

//  Setup trajectory for TL and AD
    oops::instantiateLinearModelFactory<MODEL>();
    oops::PostProcessor<State_> post;
    oops::PostProcessorTLAD<MODEL> pptraj;
    tlm_.reset(new LinearModel_(*resol_, *tlConf_));
    post.enrollProcessor(new oops::TrajectorySaver<MODEL>(*tlConf_, *resol_, *bias_, tlm_, pptraj));
    State_ xx((*xref_)[0]);
    model_->forecast(xx, *bias_, len, post);
  }

  ~LinearModelFixture<MODEL>() {}

  std::unique_ptr<const eckit::LocalConfiguration>   test_;
  std::unique_ptr<const eckit::LocalConfiguration>   tlConf_;
  std::unique_ptr<const Geometry_>       resol_;
  std::unique_ptr<const util::DateTime>  time_;
  std::unique_ptr<const oops::Variables> ctlvars_;
  std::unique_ptr<const State4D_>        xref_;
  std::unique_ptr<const Model_>          model_;
  std::unique_ptr<const ModelAux_>       bias_;
  std::unique_ptr<const ModelAuxIncr_>   dbias_;
  std::unique_ptr<const Covariance_>     B_;
  std::shared_ptr<LinearModel_>          tlm_;
};

// =============================================================================

/// \brief tests constructor, timeResolution() method and print method
template <typename MODEL> void testLinearModelConstructor() {
  typedef LinearModelFixture<MODEL>   Test_;

  const util::Duration zero(0);
  EXPECT(Test_::tlm().timeResolution() > zero);
  oops::Log::test() << "Testing LinearModel: " << Test_::tlm() << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearModelZeroLength() {
  typedef LinearModelFixture<MODEL>      Test_;
  typedef oops::Increment<MODEL>         Increment_;
  typedef oops::Increment4D<MODEL>       Increment4D_;
  typedef oops::ModelAuxIncrement<MODEL> ModelAuxIncr_;

  const util::DateTime vt(Test_::time());
  const util::Duration zero(0);

  Increment4D_ dxref(Test_::resol(), Test_::ctlvars(), {vt});
  Test_::covariance().randomize(dxref);
  ModelAuxIncr_ daux(Test_::dbias());
  const double ininorm = dxref[0].norm();
  EXPECT(ininorm > 0.0);

  Increment_ dx(dxref[0]);
  Test_::tlm().forecastTL(dx, daux, zero);
  EXPECT(dx.validTime() == vt);
  // In principle, test that a zero-length forecast (=> calling initialize then finalize) gives
  // back the initial value of the Increment. In reality, for some models, the value is changed,
  // for example if some accuracy-degrading variable change is used in going from the DA variables
  // to/from the model internal variables. To handle this case, add a tolerance:
  const double tol = Test_::test().getDouble("tolerance zero length forecast", 0.0);
  EXPECT(oops::is_close(dx.norm(), ininorm, tol));

  dx = dxref[0];
  Test_::tlm().forecastAD(dx, daux, zero);
  EXPECT(dx.validTime() == vt);
  EXPECT(oops::is_close(dx.norm(), ininorm, tol));
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearModelZeroPert() {
  typedef LinearModelFixture<MODEL>      Test_;
  typedef oops::Increment<MODEL>         Increment_;
  typedef oops::ModelAuxIncrement<MODEL> ModelAuxIncr_;

  const util::Duration len(Test_::test().getString("forecast length"));
  const util::DateTime t1(Test_::time());
  const util::DateTime t2(t1 + len);
  EXPECT(t2 > t1);

  Increment_ dx(Test_::resol(), Test_::ctlvars(), t1);
  ModelAuxIncr_ daux(Test_::dbias());

  dx.zero();
  daux.zero();
  EXPECT(dx.norm() == 0.0);
  Test_::tlm().forecastTL(dx, daux, len);
  EXPECT(dx.validTime() == t2);
  EXPECT(dx.norm() == 0.0);

  dx.zero();
  daux.zero();
  EXPECT(dx.norm() == 0.0);
  Test_::tlm().forecastAD(dx, daux, len);
  EXPECT(dx.validTime() == t1);
  EXPECT(dx.norm() == 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearModelLinearity() {
  typedef LinearModelFixture<MODEL>      Test_;
  typedef oops::Increment<MODEL>         Increment_;
  typedef oops::Increment4D<MODEL>       Increment4D_;
  typedef oops::ModelAuxIncrement<MODEL> ModelAuxIncr_;

  const util::Duration len(Test_::test().getString("forecast length"));
  const util::DateTime t1(Test_::time());
  const util::DateTime t2(t1 + len);
  EXPECT(t2 > t1);
  const double zz = 3.1415;

  Increment4D_ dx1(Test_::resol(), Test_::ctlvars(), {t1});
  Test_::covariance().randomize(dx1);
  ModelAuxIncr_ daux1(Test_::dbias());
  EXPECT(dx1[0].norm() > 0.0);

  Increment_ dx2(dx1[0]);
  ModelAuxIncr_ daux2(daux1);

  Test_::tlm().forecastTL(dx1[0], daux1, len);
  EXPECT(dx1[0].validTime() == t2);
  dx1 *= zz;
  daux1 *= zz;

  dx2 *= zz;
  daux2 *= zz;
  Test_::tlm().forecastTL(dx2, daux2, len);
  EXPECT(dx2.validTime() == t2);

  const double tol = Test_::test().getDouble("tolerance AD");
  EXPECT(oops::is_close(dx1[0].norm(), dx2.norm(), tol));
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearApproximation() {
  typedef LinearModelFixture<MODEL>      Test_;
  typedef oops::Increment<MODEL>         Increment_;
  typedef oops::Increment4D<MODEL>       Increment4D_;
  typedef oops::State<MODEL>             State_;

  const util::Duration len(Test_::test().getString("forecast length"));
  const util::DateTime t1(Test_::time());
  const util::DateTime t2(t1 + len);
  EXPECT(t2 > t1);

  Increment4D_ dx0(Test_::resol(), Test_::ctlvars(), {t1});
  Test_::covariance().randomize(dx0);
  EXPECT(dx0[0].norm() > 0.0);

  Increment_ dx(dx0[0]);
  Test_::tlm().forecastTL(dx, Test_::dbias(), len);
  const double dxnorm = dx.norm();

  oops::PostProcessor<State_> post;
  State_ xx0(Test_::xref()[0]);
  Test_::model().forecast(xx0, Test_::bias(), len, post);

  const unsigned int ntest = Test_::test().getInt("iterations TL");
  double zz = 1.0;
  if (Test_::test().has("first multiplier TL")) {
    zz = Test_::test().getDouble("first multiplier TL");
  }

  std::vector<double> errors;
  for (unsigned int jtest = 0; jtest < ntest; ++jtest) {
    State_ xx(Test_::xref()[0]);
    Increment_ pert(dx0[0]);
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
    std::streamsize ss = oops::Log::test().precision();
    oops::Log::test() << std::setprecision(16);
    oops::Log::test() << "TL error = " << err
                      << ", relative error = " << errnorm / difnorm << std::endl;
    oops::Log::test() << std::setprecision(ss);
    zz /= 10.0;
  }

// Analyze results
  const double approx = *std::min_element(errors.begin(), errors.end());
  oops::Log::test() << "Test TL min error = " << approx << std::endl;
  const double tol = Test_::test().getDouble("tolerance TL");
  EXPECT(approx < tol);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearModelAdjoint() {
  typedef LinearModelFixture<MODEL>      Test_;
  typedef oops::Increment<MODEL>         Increment_;
  typedef oops::Increment4D<MODEL>       Increment4D_;
  typedef oops::ModelAuxIncrement<MODEL> ModelAuxIncr_;

  const util::Duration len(Test_::test().getString("forecast length"));
  const util::DateTime t1(Test_::time());
  const util::DateTime t2(t1 + len);
  EXPECT(t2 > t1);

  Increment4D_ dx11(Test_::resol(), Test_::ctlvars(), {t1});
  Test_::covariance().randomize(dx11);
  ModelAuxIncr_ daux1(Test_::dbias());
  EXPECT(dx11[0].norm() > 0.0);
  Increment_ dx12(dx11[0]);
  Test_::tlm().forecastTL(dx12, daux1, len);
  EXPECT(dx12.norm() > 0.0);
  Increment4D_ dx22(Test_::resol(), Test_::ctlvars(), {t2});
  Test_::covariance().randomize(dx22);
  ModelAuxIncr_ daux2(Test_::dbias());
  EXPECT(dx22[0].norm() > 0.0);
  Increment_ dx21(dx22[0]);
  Test_::tlm().forecastAD(dx21, daux2, len);
  EXPECT(dx21.norm() > 0.0);

  EXPECT(dx11[0].norm() != dx22[0].norm());
  EXPECT(dx11[0].validTime() == t1);
  EXPECT(dx21.validTime() == t1);
  EXPECT(dx12.validTime() == t2);
  EXPECT(dx22[0].validTime() == t2);

  const double dot1 = dot_product(dx11[0], dx21);
  const double dot2 = dot_product(dx12, dx22[0]);
  const double tol = Test_::test().getDouble("tolerance AD");
  EXPECT(oops::is_close(dot1, dot2, tol));
}

// =============================================================================

template <typename MODEL>
class LinearModel : public oops::Test {
 public:
  using oops::Test::Test;
  virtual ~LinearModel() {LinearModelFixture<MODEL>::reset();}

 private:
  std::string testid() const override {return "test::LinearModel<" + MODEL::name() + ">";}

  // Override the base-class execute method to set up the LinearModelFixture by
  // calling getInstance() with the MPI commuicator required on the first call.
  // This requires that the TestEnvironment is set up beforehand.
  // Then, continue with the base-class method.
  int execute(const eckit::Configuration & globalConf, bool validate) const override {
    typedef LinearModelFixture<MODEL> Test_;

    TestEnvironment::getInstance().setup(globalConf);

    // Initialise the LinearModelFixture singleton with the communicator.
    Test_::getInstance(&getComm());

    return oops::Test::execute(globalConf, validate);
  }

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/LinearModel/testLinearModelConstructor")
      { testLinearModelConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/LinearModel/testLinearModelZeroLength")
      { testLinearModelZeroLength<MODEL>(); });
    ts.emplace_back(CASE("interface/LinearModel/testLinearModelZeroPert")
      { testLinearModelZeroPert<MODEL>(); });
    ts.emplace_back(CASE("interface/LinearModel/testLinearModelLinearity")
      { testLinearModelLinearity<MODEL>(); });
    ts.emplace_back(CASE("interface/LinearModel/testLinearApproximation")
      { testLinearApproximation<MODEL>(); });
    ts.emplace_back(CASE("interface/LinearModel/testLinearModelAdjoint")
      { testLinearModelAdjoint<MODEL>(); });
  }

  void clear() const override {}
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_LINEARMODEL_H_
