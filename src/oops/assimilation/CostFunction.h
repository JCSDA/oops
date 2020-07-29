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
#include <vector>
#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostJbState.h"
#include "oops/assimilation/CostJbTotal.h"
#include "oops/assimilation/CostJo.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/assimilation/DualVector.h"
#include "oops/assimilation/JqTerm.h"
#include "oops/assimilation/JqTermTLAD.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/TrajectorySaver.h"
#include "oops/base/VariableChangeBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/LinearModel.h"
#include "oops/interface/State.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

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
  typedef JqTerm<MODEL>                 JqTerm_;
  typedef JqTermTLAD<MODEL>             JqTermTLAD_;
  typedef Geometry<MODEL>               Geometry_;
  typedef LinearModel<MODEL>            LinearModel_;
  typedef State<MODEL>                  State_;
  typedef Increment<MODEL>              Increment_;

 public:
  explicit CostFunction(const eckit::Configuration &);
  virtual ~CostFunction() {}

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
  const CostBase_ & jterm(const unsigned ii) const {return jterms_[ii];}
  unsigned nterms() const {return jterms_.size();}
  double getCostJb() const {return costJb_;}
  double getCostJoJc() const {return costJoJc_;}

 protected:
  void setupTerms(const eckit::Configuration &);
  void setupTerms(const eckit::Configuration &, const State_ &);
  const LinearModel_ & getTLM(const unsigned isub = 0) const {return tlm_[isub];}
  const CtrlVar_ & background() const {return *xb_;}

 private:
  virtual void addIncr(CtrlVar_ &, const CtrlInc_ &, PostProcessor<Increment_>&) const = 0;

  virtual CostJbState<MODEL>  * newJb(const eckit::Configuration &, const Geometry_ &,
                                      const CtrlVar_ &) const = 0;
  virtual CostJo<MODEL, OBS>       * newJo(const eckit::Configuration &) const = 0;
  virtual CostTermBase<MODEL, OBS> * newJc(const eckit::Configuration &,
                                           const Geometry_ &) const = 0;
  virtual void doLinearize(const Geometry_ &, const eckit::Configuration &,
                           const CtrlVar_ &, const CtrlVar_ &) = 0;
  virtual const Geometry_ & geometry() const = 0;

// Data members
  std::unique_ptr<const CtrlVar_> xb_;
  std::unique_ptr<JbTotal_> jb_;
  boost::ptr_vector<CostBase_> jterms_;
  boost::ptr_vector<LinearModel_> tlm_;

  mutable double costJb_;
  mutable double costJoJc_;
};

// -----------------------------------------------------------------------------

/// Cost Function Factory
template <typename MODEL, typename OBS>
class CostFactory {
 public:
  static CostFunction<MODEL, OBS> * create(const eckit::Configuration &, const eckit::mpi::Comm &);
  virtual ~CostFactory() = default;

 protected:
  explicit CostFactory(const std::string &);
 private:
  virtual CostFunction<MODEL, OBS> * make(const eckit::Configuration &,
                                          const eckit::mpi::Comm &) = 0;
  static std::map < std::string, CostFactory<MODEL, OBS> * > & getMakers() {
    static std::map < std::string, CostFactory<MODEL, OBS> * > makers_;
    return makers_;
  }
};

template<class MODEL, class OBS, class FCT>
class CostMaker : public CostFactory<MODEL, OBS> {
 private:
  CostFunction<MODEL, OBS> * make(const eckit::Configuration & config,
                             const eckit::mpi::Comm & comm) override
    {return new FCT(config, comm);}
 public:
  explicit CostMaker(const std::string & name) : CostFactory<MODEL, OBS>(name) {}
};

// =============================================================================

//  Factory
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostFactory<MODEL, OBS>::CostFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    Log::error() << name << " already registered in cost function factory." << std::endl;
    ABORT("Element already registered in CostFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostFunction<MODEL, OBS>* CostFactory<MODEL, OBS>::create(const eckit::Configuration & config,
                                                const eckit::mpi::Comm & comm) {
  std::string id = config.getString("cost type");
  Log::trace() << "Variational Assimilation Type=" << id << std::endl;
  typename std::map<std::string, CostFactory<MODEL, OBS>*>::iterator j = getMakers().find(id);
  if (j == getMakers().end()) {
    Log::error() << id << " does not exist in cost function factory." << std::endl;
    ABORT("Element does not exist in CostFactory.");
  }
  Log::trace() << "CostFactory::create found cost function type" << std::endl;
  return (*j).second->make(config, comm);
}

// -----------------------------------------------------------------------------
//  Cost Function
// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
CostFunction<MODEL, OBS>::CostFunction(const eckit::Configuration & config)
  : jb_(), jterms_(), tlm_()
{}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFunction<MODEL, OBS>::setupTerms(const eckit::Configuration & config) {
  Log::trace() << "CostFunction::setupTerms start" << std::endl;

// Jo
  CostJo<MODEL, OBS> * jo = this->newJo(config);
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
void CostFunction<MODEL, OBS>::setupTerms(const eckit::Configuration & config,
                                          const State_ & statein) {
  Log::trace() << "CostFunction::setupTerms start" << std::endl;

// Jo
  CostJo<MODEL, OBS> * jo = this->newJo(config);
  jterms_.push_back(jo);
  Log::trace() << "CostFunction::setupTerms Jo added" << std::endl;

// Jb
  xb_.reset(new CtrlVar_(config, statein, jo->obspaces()));
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
  PostProcessor<State_> pp(post);
  JqTerm_ * jq = jb_->initialize(fguess);
  pp.enrollProcessor(jq);
  for (unsigned jj = 0; jj < jterms_.size(); ++jj) {
    pp.enrollProcessor(jterms_[jj].initialize(fguess, config));
  }

// Run NL model
  CtrlVar_ mfguess(fguess);
  this->runNL(mfguess, pp);

// Cost function value
  double zzz = 0.0;
  costJb_ = jb_->finalize(jq);
  zzz += costJb_;
  costJoJc_ = 0.0;
  for (unsigned jj = 0; jj < jterms_.size(); ++jj) {
    costJoJc_ += jterms_[jj].finalize();
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
  const Geometry_ lowres(resConf, this->geometry().getComm());

// Setup trajectory for terms of cost function
  PostProcessorTLAD<MODEL> pptraj;
  JqTermTLAD_ * jq = jb_->initializeTraj(fguess, lowres);
  pptraj.enrollProcessor(jq);
  for (unsigned jj = 0; jj < jterms_.size(); ++jj) {
    pptraj.enrollProcessor(jterms_[jj].initializeTraj(fguess, lowres, innerConf));
  }

// Setup linear model (and trajectory)
// YT: TrajectorySaver should be QuadraticCostFunction
  const eckit::LocalConfiguration tlmConf(innerConf, "linear model");
  tlm_.clear();   // YT: Should release at the end and should be inside quadratic J object
  PostProcessor<State_> pp(post);
  pp.enrollProcessor(new TrajectorySaver<MODEL>(tlmConf, lowres, fguess.modVar(), tlm_, pptraj));

// Run NL model
  double zzz = this->evaluate(fguess, innerConf, pp);

// Finalize trajectory setup
  jb_->finalizeTraj(jq);
  for (unsigned jj = 0; jj < jterms_.size(); ++jj) {
    jterms_[jj].finalizeTraj();
  }

// Specific linearization if needed
  this->doLinearize(lowres, innerConf, *xb_, fguess);

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

  for (unsigned jj = 0; jj < jterms_.size(); ++jj) {
    std::shared_ptr<const GeneralizedDepartures> tmp(jterms_[jj].newGradientFG());
    costad.enrollProcessor(jterms_[jj].setupAD(tmp, grad));
  }

  this->runADJ(grad, costad, pp);
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
  for (unsigned jj = 0; jj < jterms_.size(); ++jj) {
    jterms_[jj].resetLinearization();
  }
  tlm_.clear();
  Log::trace() << "CostFunction::resetLinearization done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTFUNCTION_H_
