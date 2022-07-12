/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_COSTFUNCTION_H_
#define OOPS_ASSIMILATION_COSTFUNCTION_H_

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostJbTotal.h"
#include "oops/assimilation/CostJo.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/ObsTypeParameters.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"

namespace oops {

// Forward decl
template <typename MODEL, typename OBS> class CostFactory;
template <typename MODEL> class CostJbState;

// -----------------------------------------------------------------------------

/// \brief Parameters for an implementation of CostFunction
template <typename MODEL, typename OBS>
class CostFunctionParametersBase : public Parameters {
  OOPS_ABSTRACT_PARAMETERS(CostFunctionParametersBase, Parameters)

  typedef typename Geometry<MODEL>::Parameters_ GeometryParameters_;

 public:
  RequiredParameter<std::string> costType{"cost type", "Name of the cost function", this};

  // These options aren't used in the CostFunction base class, but are used in the constructors of
  // each individual cost function implementation.
  RequiredParameter<GeometryParameters_> geometry{"geometry", "geometry", this};
  RequiredParameter<util::Duration> windowLength{"window length", "length of assimilation window",
      this};
  RequiredParameter<util::DateTime> windowBegin{"window begin",
      "start datetime of assimilation window", this};
  RequiredParameter<Variables> analysisVariables{"analysis variables", "variables to assimilate",
      this};

  // These options are used (directly or indirectly) in the CostFunction base class.
  // Jo
  RequiredParameter<ObserversParameters<MODEL, OBS>> observations{"observations", this};

  // Jb
  Parameter<bool> jbEvaluation{"jb evaluation", true, this};

  // Jc
  OptionalParameter<eckit::LocalConfiguration> constraints{"constraints", this};

  // Some models use extra options, which we will eventually clean up. Until then, make tests pass
  // with this config parameter to "catch" options not described by Parameters.
  ConfigurationParameter extraOptionsInConfig{this};
};

// -----------------------------------------------------------------------------

/// \brief Contains a polymorphic parameter holding an instance of a subclass of
/// CostFunctionParametersBase.
template <typename MODEL, typename OBS>
class CostFunctionParametersWrapper : public Parameters {
  OOPS_CONCRETE_PARAMETERS(CostFunctionParametersWrapper, Parameters)

 public:
  /// After deserialization, holds an instance of a subclass of CostFunctionParametersBase
  /// controlling the behavior of a cost function. The type of the subclass is
  /// determined by the value of the "cost type" key in the Configuration object from which
  /// this object is deserialized.
  RequiredPolymorphicParameter<CostFunctionParametersBase<MODEL, OBS>, CostFactory<MODEL, OBS>>
    costTypeParameters{"cost type", "Name of the cost function", this};
};

// -----------------------------------------------------------------------------

/// Cost Function
/*!
 * The CostFunction defines and manages the computation of all the terms
 * of the variational data assimilation cost function.
 */
template<typename MODEL, typename OBS> class CostFunction : private boost::noncopyable {
  typedef ControlIncrement<MODEL, OBS>  CtrlInc_;
  typedef ControlVariable<MODEL, OBS>   CtrlVar_;
  typedef CostJbTotal<MODEL, OBS>       JbTotal_;
  typedef CostTermBase<MODEL, OBS>      CostBase_;
  typedef Geometry<MODEL>               Geometry_;
  typedef State<MODEL>                  State_;
  typedef Increment<MODEL>              Increment_;

 public:
  typedef CostFunctionParametersWrapper<MODEL, OBS> Parameters_;

  CostFunction() = default;
  virtual ~CostFunction() = default;

  double evaluate(const CtrlVar_ &,
                  const eckit::Configuration & config = eckit::LocalConfiguration(),
                  PostProcessor<State_> post = PostProcessor<State_>() );
  double linearize(const CtrlVar_ &, const eckit::Configuration &,
                   PostProcessor<State_> post = PostProcessor<State_>() );

  virtual void runTLM(CtrlInc_ &, PostProcessorTLAD<MODEL> &,
                      PostProcessor<Increment_> post = PostProcessor<Increment_>(),
                      const bool idModel = false) const = 0;
  virtual void runADJ(CtrlInc_ &, PostProcessorTLAD<MODEL> &,
                      PostProcessor<Increment_> post = PostProcessor<Increment_>(),
                      const bool idModel = false) const = 0;
  virtual void zeroAD(CtrlInc_ &) const = 0;

  virtual void runNL(CtrlVar_ &, PostProcessor<State_>&) const = 0;

  void addIncrement(CtrlVar_ &, const CtrlInc_ &,
                    PostProcessor<Increment_> post = PostProcessor<Increment_>() ) const;
  void resetLinearization();

/// Compute cost function gradient at first guess (without Jb).
  void computeGradientFG(CtrlInc_ &) const;

/// Access \f$ J_b\f$
  const JbTotal_ & jb() const {return *jb_;}
/// Access terms of the cost function other than \f$ J_b\f$
  const CostBase_ & jterm(const size_t ii) const {return jterms_[ii];}
  size_t nterms() const {return jterms_.size();}
  double getCostJb() const {return costJb_;}
  double getCostJoJc() const {return costJoJc_;}

 protected:
  void setupTerms(const eckit::Configuration &);
  const CtrlVar_ & background() const {return *xb_;}

 private:
  virtual void addIncr(CtrlVar_ &, const CtrlInc_ &, PostProcessor<Increment_>&) const = 0;

  virtual CostJbState<MODEL>  * newJb(const eckit::Configuration &, const Geometry_ &,
                                      const CtrlVar_ &) const = 0;
  virtual CostJo<MODEL, OBS>       * newJo(const ObserversParameters<MODEL, OBS> &) const = 0;
  virtual CostTermBase<MODEL, OBS> * newJc(const eckit::Configuration &,
                                           const Geometry_ &) const = 0;
  virtual void doLinearize(const Geometry_ &, const eckit::Configuration &,
                           const CtrlVar_ &, const CtrlVar_ &,
                           PostProcessor<State_> &, PostProcessorTLAD<MODEL> &) = 0;
  virtual const Geometry_ & geometry() const = 0;

// Data members
  std::unique_ptr<const CtrlVar_> xb_;
  std::unique_ptr<JbTotal_> jb_;
  boost::ptr_vector<CostBase_> jterms_;
  std::unique_ptr<const Geometry_> lowres_;
  // Geometry from the previous iteration. A temporary fix Required because
  // lowres_ is not properly shared between all the objects that use it.
  // lowres_ gets destructed after the first outer loop iteration, but some
  // minimizers need access to the first outer loop Geometry in the second
  // outer loop (e.g. DRMinimizer).
  std::unique_ptr<const Geometry_> lowres_previter_;

  mutable double costJb_;
  mutable double costJoJc_;
};

// -----------------------------------------------------------------------------

/// Cost Function Factory
template <typename MODEL, typename OBS>
class CostFactory {
 public:
  static std::unique_ptr<CostFunction<MODEL, OBS>> create(
      const CostFunctionParametersBase<MODEL, OBS> &, const eckit::mpi::Comm &);
  virtual ~CostFactory() = default;

  static std::unique_ptr<CostFunctionParametersBase<MODEL, OBS>> createParameters(
      const std::string &name);

  static std::vector<std::string> getMakerNames() {
    return keys(getMakers());
  }

 protected:
  explicit CostFactory(const std::string &);

 private:
  virtual std::unique_ptr<CostFunction<MODEL, OBS>> make(
      const CostFunctionParametersBase<MODEL, OBS> & params, const eckit::mpi::Comm &) = 0;

  virtual std::unique_ptr<CostFunctionParametersBase<MODEL, OBS>> makeParameters() const = 0;

  static std::map < std::string, CostFactory<MODEL, OBS> * > & getMakers() {
    static std::map < std::string, CostFactory<MODEL, OBS> * > makers_;
    return makers_;
  }
};

template<class MODEL, class OBS, class FCT>
class CostMaker : public CostFactory<MODEL, OBS> {
 private:
  typedef typename FCT::Parameters_ Parameters_;
  std::unique_ptr<CostFunction<MODEL, OBS>> make(
      const CostFunctionParametersBase<MODEL, OBS> & params,
      const eckit::mpi::Comm & comm) override
  {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return std::make_unique<FCT>(stronglyTypedParams, comm);
  }

  std::unique_ptr<CostFunctionParametersBase<MODEL, OBS>> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit CostMaker(const std::string & name) : CostFactory<MODEL, OBS>(name) {}
};

// =============================================================================

//  Factory
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostFactory<MODEL, OBS>::CostFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in cost function factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
std::unique_ptr<CostFunctionParametersBase<MODEL, OBS>>
CostFactory<MODEL, OBS>::createParameters(const std::string &name) {
  typename std::map<std::string, CostFactory<MODEL, OBS>*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in CostFactory");
  }
  return it->second->makeParameters();
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
std::unique_ptr<CostFunction<MODEL, OBS>> CostFactory<MODEL, OBS>::create(
    const CostFunctionParametersBase<MODEL, OBS> & params, const eckit::mpi::Comm & comm) {
  const std::string id = params.costType;
  Log::trace() << "Variational Assimilation Type=" << id << std::endl;
  typename std::map<std::string, CostFactory<MODEL, OBS>*>::iterator j = getMakers().find(id);
  if (j == getMakers().end()) {
    throw std::runtime_error(id + " does not exist in cost function factory.");
  }
  std::unique_ptr<CostFunction<MODEL, OBS>> ptr(j->second->make(params, comm));
  Log::trace() << "CostFactory::create found cost function type" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------
//  Cost Function
// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFunction<MODEL, OBS>::setupTerms(const eckit::Configuration & config) {
  Log::trace() << "CostFunction::setupTerms start" << std::endl;

// Jo
  eckit::LocalConfiguration obsconf(config, "observations");
  ObserversParameters<MODEL, OBS> joparams{};
  joparams.deserialize(obsconf);
  CostJo<MODEL, OBS> * jo = this->newJo(joparams);
  jterms_.push_back(jo);
  Log::trace() << "CostFunction::setupTerms Jo added" << std::endl;

// Jb
  xb_.reset(new CtrlVar_(config, this->geometry(), jo->obspaces()));
  jb_.reset(new JbTotal_(*xb_, this->newJb(config, this->geometry(), *xb_),
                         config, this->geometry(), jo->obspaces()));
  Log::trace() << "CostFunction::setupTerms Jb added" << std::endl;

// Other constraints
  std::vector<eckit::LocalConfiguration> jcs;
  config.get("constraints", jcs);
  for (size_t jj = 0; jj < jcs.size(); ++jj) {
    CostTermBase<MODEL, OBS> * jc = this->newJc(jcs[jj], this->geometry());
    jterms_.push_back(jc);
  }
  Log::trace() << "CostFunction::setupTerms done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
double CostFunction<MODEL, OBS>::evaluate(const CtrlVar_ & fguess,
                                          const eckit::Configuration & config,
                                          PostProcessor<State_> post) {
  Log::trace() << "CostFunction::evaluate start" << std::endl;
// Setup terms of cost function
  jb_->initialize(fguess);
  PostProcessor<State_> pp(post);
  for (size_t jj = 0; jj < jterms_.size(); ++jj) {
    jterms_[jj].setPostProc(fguess, config, pp);
  }

// Run NL model
  CtrlVar_ mfguess(fguess);
  this->runNL(mfguess, pp);

// Cost function value
  double zzz = 0.0;
  costJb_ = jb_->finalize(mfguess);
  zzz += costJb_;
  costJoJc_ = 0.0;
  for (size_t jj = 0; jj < jterms_.size(); ++jj) {
    costJoJc_ += jterms_[jj].computeCost();
  }
  zzz += costJoJc_;
  Log::test() << "CostFunction: Nonlinear J = " << zzz << std::endl;
  Log::trace() << "CostFunction::evaluate done" << std::endl;
  return zzz;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
double CostFunction<MODEL, OBS>::linearize(const CtrlVar_ & fguess,
                                           const eckit::Configuration & innerConf,
                                           PostProcessor<State_> post) {
  Log::trace() << "CostFunction::linearize start" << std::endl;
// Inner loop resolution
  const eckit::LocalConfiguration resConf(innerConf, "geometry");
  lowres_previter_ = std::move(lowres_);
  lowres_ = std::make_unique<Geometry_>(resConf, this->geometry().getComm(),
                                        this->geometry().timeComm());

// Setup trajectory for terms of cost function
  PostProcessorTLAD<MODEL> pptraj;
  jb_->initializeTraj(fguess, *lowres_, innerConf, pptraj);
  for (size_t jj = 0; jj < jterms_.size(); ++jj) {
    jterms_[jj].setPostProcTraj(fguess, innerConf, *lowres_, pptraj);
  }

// Specific linearization if needed (including TLM)
  this->doLinearize(*lowres_, innerConf, *xb_, fguess, post, pptraj);

// Run NL model
  double zzz = this->evaluate(fguess, innerConf, post);

// Finalize trajectory setup
  jb_->finalizeTraj();
  for (size_t jj = 0; jj < jterms_.size(); ++jj) {
    jterms_[jj].computeCostTraj();
  }

  Log::trace() << "CostFunction::linearize done" << std::endl;
  return zzz;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFunction<MODEL, OBS>::computeGradientFG(CtrlInc_ & grad) const {
  Log::trace() << "CostFunction::computeGradientFG start" << std::endl;
  PostProcessor<Increment_> pp;
  PostProcessorTLAD<MODEL> costad;
  this->zeroAD(grad);

  for (size_t jj = 0; jj < jterms_.size(); ++jj) {
    std::shared_ptr<const GeneralizedDepartures> tmp(jterms_[jj].newGradientFG());
    jterms_[jj].computeCostAD(tmp, grad, costad);
  }

  this->runADJ(grad, costad, pp);

  for (size_t jj = 0; jj < jterms_.size(); ++jj) {
    jterms_[jj].setPostProcAD();
  }

  Log::info() << "CostFunction::computeGradientFG: gradient:" << grad << std::endl;
  Log::trace() << "CostFunction::computeGradientFG done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFunction<MODEL, OBS>::addIncrement(CtrlVar_ & xx, const CtrlInc_ & dx,
                                            PostProcessor<Increment_> post) const {
  Log::trace() << "CostFunction::addIncrement start" << std::endl;
  Log::info() << "CostFunction::addIncrement: First guess:" << xx << std::endl;
  Log::info() << "CostFunction::addIncrement: Increment:" << dx << std::endl;

  xx.obsVar() += dx.obsVar();
  xx.modVar() += dx.modVar();
  this->addIncr(xx, dx, post);

  Log::info() << "CostFunction::addIncrement: Analysis: " << xx << std::endl;
  Log::test() << "CostFunction::addIncrement: Analysis: " << xx << std::endl;
  Log::trace() << "CostFunction::addIncrement done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFunction<MODEL, OBS>::resetLinearization() {
  Log::trace() << "CostFunction::resetLinearization start" << std::endl;
  for (size_t jj = 0; jj < jterms_.size(); ++jj) {
    jterms_[jj].resetLinearization();
  }
  Log::trace() << "CostFunction::resetLinearization done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTFUNCTION_H_
