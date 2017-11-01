/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_BASE_TESTSUITEOPOBSFIXTURE_H_
#define TEST_BASE_TESTSUITEOPOBSFIXTURE_H_

#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/pointer_cast.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/test/unit_test.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "util/Logger.h"
#include "util/DateTime.h"
#include "util/Duration.h"

#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/JqTerm.h"
#include "oops/assimilation/DualVector.h"
#include "oops/assimilation/HMatrix.h"
#include "oops/assimilation/HtMatrix.h"
#include "oops/base/Departures.h"
#include "oops/base/Observations.h"
#include "oops/base/Observer.h"
#include "oops/base/ObserverTL.h"
#include "oops/base/ObserverAD.h"
#include "oops/base/ModelIncrement.h"
#include "oops/base/ModelState.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTL.h"
#include "oops/base/PostProcessorAD.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/GeneralizedDepartures.h"
#include "oops/base/StateInfo.h"
#include "oops/interface/ObsErrorCovariance.h"

#include "oops/assimilation/instantiateCostFactory.h"
#include "oops/assimilation/instantiateMinFactory.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/instantiateObsErrorFactory.h"

namespace test {

/*!
 *  This class gives access to input datas for the whole test suite. These
 *  datas are calculated once.
 */
template<typename MODEL> class TestSuiteOpObsFixture {
  public:

    typedef typename MODEL::Geometry Geometry_;
    typedef typename MODEL::ModelAuxControl ModelAuxCtrl_;
    typedef typename MODEL::ModelAuxCovariance ModelAuxCovar_;
    typedef typename MODEL::ModelAuxIncrement ModelAuxIncrement_;
    typedef typename MODEL::ModelConfiguration ModelConfig_;
    typedef typename MODEL::LinearModel LinearModel_;
    typedef typename MODEL::ObsAuxIncrement ObsAuxIncr_;
    typedef typename MODEL::ObsAuxControl ObsAuxCtrl_;
    typedef typename MODEL::ObsAuxCovariance ObsAuxCovar_;
    typedef typename MODEL::State State_;
    typedef typename MODEL::Increment Increment_;
    typedef typename MODEL::ObsOperatorLinearizationTrajectory ObsOpTraj_;
    typedef typename MODEL::ModelAuxIncrement ModelAuxIncr_;
    typedef oops::ModelState<MODEL> ModelState_;
    typedef oops::Departures<MODEL> Departures_;
    typedef oops::Observations<MODEL> Observations_;
    typedef oops::ObsErrorCovariance<MODEL> ObsErrorCovariance_;
    typedef oops::PostProcessor<State_> PostProcessor_;
    typedef oops::PostProcessorTL<Increment_> PostProcessorTL_;
    typedef oops::PostProcessorAD<Increment_> PostProcessorAD_;
    typedef oops::ObserverTL<MODEL, Increment_> ObserverTL_;
    typedef oops::ModelIncrement<MODEL> ModelIncrement_;
    typedef oops::ModelSpaceCovarianceBase<MODEL> ModelSpaceCovarianceBase_;
    typedef oops::ControlIncrement<MODEL> CtrlInc_;
    typedef oops::ControlVariable<MODEL> ControlVariable_;
    typedef oops::DualVector<MODEL> DepartVector_;
    typedef oops::HMatrix<MODEL> HMatrix_;
    typedef oops::HtMatrix<MODEL> HtMatrix_;
    typedef oops::CostFunction<MODEL> CostFunction_;
    typedef oops::CostFactory<MODEL> CostFactory_;
    typedef oops::Observer<MODEL, State_> Observer_;

    // -----------------------------------------------------------------------------
    TestSuiteOpObsFixture() : tolAD_(0.0), deltaX_(1.0) {}
    // -----------------------------------------------------------------------------
    ~TestSuiteOpObsFixture() {}
    // -----------------------------------------------------------------------------
    const CostFunction_ & J() {
      return *getInstance().J_;
    }

    // -----------------------------------------------------------------------------
    double tolAD() {
      return getInstance().tolAD_;
    }

    // -----------------------------------------------------------------------------
    double deltaX() {
      return getInstance().deltaX_;
    }

    // -----------------------------------------------------------------------------
    const util::Config & getConfig() {
      return *getInstance().config_;
    }

    // -----------------------------------------------------------------------------
    const ModelConfig_ & getModelConfig() {
      return *getInstance().model_;
    }

    // -----------------------------------------------------------------------------
    const ModelAuxCtrl_ & getModelErr() {
      return *getInstance().moderr_;
    }

    // -----------------------------------------------------------------------------
    boost::shared_ptr<Departures_> getRandomDepartures(unsigned jj) {
      // Cannot randomize dy directly for matR is hidden within a costJo object
      boost::scoped_ptr<CtrlInc_> dx(new CtrlInc_(J().jb()));
      randomize(*dx);

      PostProcessorTL_ cost;
      cost.enrollProcessor(J().jterm(jj).setupTL(*dx));

      J().runTLM(*dx, cost);
      boost::shared_ptr<Departures_> dy(
          boost::dynamic_pointer_cast<Departures_>(
              cost.releaseOutputFromTL(0)));

      (*dy) *= 3.0;
      return dy;
    }

    // -----------------------------------------------------------------------------
    void randomize(DepartVector_ & dy) {
      // Cannot randomize dy directly for matR is hidden within a costJo object
      const HMatrix_ H(J());

      boost::scoped_ptr<CtrlInc_> dx(new CtrlInc_(J().jb()));
      randomize(*dx);

      dy.clear();
      H.multiply(*dx, dy);

      dy *= 3.0;
    }

    // -----------------------------------------------------------------------------
    void randomize(CtrlInc_ & dx) {
      getInstance().matB_->randomize(dx.state()[0]);
    }

    // -----------------------------------------------------------------------------
    boost::shared_ptr<Observations_> multiply_Hnl(const ControlVariable_ & xx,
        unsigned jj) {
      oops::ModelState<MODEL> xstart(xx.state()[0], *getInstance().model_);
      PostProcessor_ post;
      boost::shared_ptr<Observer_> obs(
          boost::dynamic_pointer_cast<Observer_>(J().jterm(jj).initialize(xx)));
      post.enrollProcessor(obs);

      xstart.forecast(*getInstance().moderr_, *getInstance().fclen_, post);

      boost::shared_ptr<Observations_> yequ(obs->release());

      return yequ;
    }

    // -----------------------------------------------------------------------------
    boost::shared_ptr<Departures_> multiply_Htl(const CtrlInc_ & dx,
        unsigned jj) {
      PostProcessorTL_ cost;
      cost.enrollProcessor(J().jterm(jj).setupTL(dx));

      J().runTLM(dx, cost);
      boost::shared_ptr<Departures_> dy(
          boost::dynamic_pointer_cast<Departures_>(
              cost.releaseOutputFromTL(0)));

      return dy;
    }

    // -----------------------------------------------------------------------------
    void setup(const eckit::Configuration & fullConfig) {
      config_.reset(new eckit::LocalConfiguration(fullConfig));
      tolAD_ = config_->getDouble("toleranceTestAD");
      deltaX_ = config_->getDouble("delta_x");

      oops::instantiateCostFactory<MODEL>();
      oops::instantiateCovarFactory<MODEL>();
      oops::instantiateMinFactory<MODEL>();
      oops::instantiateObsErrorFactory<MODEL>();

      //  Setup resolution
      const eckit::LocalConfiguration resolConfig(*config_, "resolution");
      const Geometry_ resol(resolConfig);

      //  Setup ModelConfig_
      const eckit::LocalConfiguration modelConfig(*config_, "model");
      model_.reset(new ModelConfig_(resol, modelConfig));

      // Setup J
      const eckit::LocalConfiguration cfConf(*config_, "cost_function");
      fclen_.reset(new util::Duration(cfConf.getString("window_length")));
      J_.reset(CostFactory_::create(cfConf, *model_));
      ControlVariable_ xx(J_->jb().getBackground());
      J_->linearize(xx, eckit::LocalConfiguration(*config_, "linear"));

      const eckit::LocalConfiguration jbConf(cfConf, "Jb");
      moderr_.reset(
          new ModelAuxCtrl_(eckit::LocalConfiguration(jbConf, "Background"), resol));

      const eckit::LocalConfiguration errCovConfig(jbConf, "Covariance");
      matB_.reset(oops::CovarianceFactory<MODEL>::create(errCovConfig, resol));
      matB_->linearize(xx.state()[0], resol);
    }

  private:

    /* This singleton must be attribute of the Run object, so that
     * all objects will be destroyed in the right order. Especially before
     * MPI finalizing. */
    TestSuiteOpObsFixture & getInstance();

    boost::scoped_ptr<util::Duration> fclen_;
    boost::scoped_ptr<eckit::LocalConfiguration> config_;
    double tolAD_, deltaX_;
    boost::scoped_ptr<ModelSpaceCovarianceBase_> matB_;
    boost::scoped_ptr<CostFunction_> J_;
    boost::scoped_ptr<ModelConfig_> model_;
    boost::scoped_ptr<ModelAuxCtrl_> moderr_;
};

template<typename MODEL>
TestSuiteOpObsFixture<MODEL> & TestSuiteOpObsFixture<MODEL>::getInstance() {
  extern TestSuiteOpObsFixture<MODEL>& getGlobalFixture();
  return getGlobalFixture();
}

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_BASE_TESTSUITEOPOBSFIXTURE_H_
