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
#include <string>
#include <vector>
#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "util/Logger.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostJbState.h"
#include "oops/assimilation/CostJbTotal.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/assimilation/CostJo.h"
#include "oops/assimilation/DualVector.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTL.h"
#include "oops/base/PostProcessorAD.h"
#include "oops/base/TrajectorySaver.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/LinearModel.h"
#include "oops/interface/Model.h"
#include "oops/interface/State.h"
#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/abor1_cpp.h"
#include "util/dot_product.h"

namespace oops {
  template<typename MODEL> class JqTerm;

// -----------------------------------------------------------------------------

/// Cost Function
/*!
 * The CostFunction defines and manages the computation of all the terms
 * of the variational data assimilation cost function.
 */

template<typename MODEL> class CostFunction : private boost::noncopyable {
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef ControlVariable<MODEL>     CtrlVar_;
  typedef CostJbTotal<MODEL>         JbTotal_;
  typedef CostTermBase<MODEL>        CostBase_;
  typedef JqTerm<MODEL>              JqTerm_;
  typedef Geometry<MODEL>            Geometry_;
  typedef Model<MODEL>               Model_;
  typedef LinearModel<MODEL>         LinearModel_;
  typedef State<MODEL>               State_;
  typedef Increment<MODEL>           Increment_;

 public:
  CostFunction(const Geometry_ &, const Model_ &);
  virtual ~CostFunction() {}

  double evaluate(const CtrlVar_ &,
                  const eckit::Configuration & config = eckit::LocalConfiguration(),
                  PostProcessor<State_> post = PostProcessor<State_>() ) const;
  double linearize(const CtrlVar_ &, const eckit::Configuration &,
                   PostProcessor<State_> post = PostProcessor<State_>() );

  virtual void runTLM(CtrlInc_ &, PostProcessorTL<Increment_> &,
                      PostProcessor<Increment_> post = PostProcessor<Increment_>(),
                      const bool idModel = false) const =0;
  virtual void runADJ(CtrlInc_ &, PostProcessorAD<Increment_> &,
                      PostProcessor<Increment_> post = PostProcessor<Increment_>(),
                      const bool idModel = false) const =0;
  virtual void zeroAD(CtrlInc_ &) const =0;

  virtual void runNL(CtrlVar_ &, PostProcessor<State_>&) const =0;

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
  const double getCostJb() const {return costJb_;}
  const double getCostJoJc() const {return costJoJc_;}

 protected:
  void setupTerms(const eckit::Configuration &);
  const Model_ & getModel() const {return model_;}
  const LinearModel_ & getTLM(const unsigned isub = 0) const {return tlm_[isub];}

 private:
  virtual void addIncr(CtrlVar_ &, const CtrlInc_ &,
                       PostProcessor<Increment_>&) const =0;

  virtual CostJbState<MODEL>  * newJb(const eckit::Configuration &, const Geometry_ &, const CtrlVar_ &) const =0;
  virtual CostJo<MODEL>       * newJo(const eckit::Configuration &) const =0;
  virtual CostTermBase<MODEL> * newJc(const eckit::Configuration &, const Geometry_ &) const =0;

// Data members
  const Geometry_ & resol_;
  const Model_ & model_;
  boost::scoped_ptr<const CtrlVar_> xb_;
  boost::scoped_ptr<JbTotal_> jb_;
  boost::ptr_vector<CostBase_> jterms_;
  boost::ptr_vector<LinearModel_> tlm_;

  double costJb_;
  double costJoJc_;
};

// -----------------------------------------------------------------------------

/// Cost Function Factory
template <typename MODEL>
class CostFactory {
  typedef Geometry<MODEL>            Geometry_;
  typedef Model<MODEL>               Model_;
 public:
  static CostFunction<MODEL> * create(const eckit::Configuration &,
                                      const Geometry_ &, const Model_ &);
  virtual ~CostFactory() { getMakers().clear(); }

 protected:
  explicit CostFactory(const std::string &);
 private:
  virtual CostFunction<MODEL> * make(const eckit::Configuration &,
                                     const Geometry_ &, const Model_ &) =0;
  static std::map < std::string, CostFactory<MODEL> * > & getMakers() {
    static std::map < std::string, CostFactory<MODEL> * > makers_;
    return makers_;
  }
};

template<class MODEL, class FCT>
class CostMaker : public CostFactory<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef Model<MODEL>               Model_;
 private:
  virtual CostFunction<MODEL> * make(const eckit::Configuration & config,
                                     const Geometry_ & resol, const Model_ & model)
    {return new FCT(config, resol, model);}
 public:
  explicit CostMaker(const std::string & name) : CostFactory<MODEL>(name) {}
};

// =============================================================================

//  Factory
// -----------------------------------------------------------------------------

template <typename MODEL>
CostFactory<MODEL>::CostFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    Log::error() << name << " already registered in cost function factory." << std::endl;
    ABORT("Element already registered in CostFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
CostFunction<MODEL>* CostFactory<MODEL>::create(const eckit::Configuration & config,
                                                const Geometry_ & resol, const Model_ & model) {
  std::string id = config.getString("cost_type");
  Log::trace() << "Variational Assimilation Type=" << id << std::endl;
  typename std::map<std::string, CostFactory<MODEL>*>::iterator j = getMakers().find(id);
  if (j == getMakers().end()) {
    Log::error() << id << " does not exist in cost function factory." << std::endl;
    ABORT("Element does not exist in CostFactory.");
  }
  return (*j).second->make(config, resol, model);
}

// -----------------------------------------------------------------------------
//  Cost Function
// -----------------------------------------------------------------------------

template<typename MODEL>
CostFunction<MODEL>::CostFunction(const Geometry_ & resol, const Model_ & model)
  : resol_(resol), model_(model), jb_(), jterms_(), tlm_()
{
  Log::trace() << "CostFunction:created" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostFunction<MODEL>::setupTerms(const eckit::Configuration & config) {
  Log::trace() << "CostFunction: setupTerms starting" << std::endl;

// Jo
  std::vector<eckit::LocalConfiguration> jos;
  config.get("Jo", jos);
  for (size_t jj = 0; jj < jos.size(); ++jj) {
    CostJo<MODEL> * jo = this->newJo(jos[jj]);
    jterms_.push_back(jo);
  }
  Log::trace() << "CostFunction: setupTerms Jo added" << std::endl;

// Jb
  const eckit::LocalConfiguration jbConf(config, "Jb");
  xb_.reset(new CtrlVar_(eckit::LocalConfiguration(jbConf, "Background"), resol_));
  jb_.reset(new JbTotal_(*xb_, this->newJb(jbConf, resol_, *xb_), jbConf, resol_));
  Log::trace() << "CostFunction: setupTerms Jb added" << std::endl;

// Other constraints
  std::vector<eckit::LocalConfiguration> jcs;
  config.get("Jc", jcs);
  for (size_t jj = 0; jj < jcs.size(); ++jj) {
    CostTermBase<MODEL> * jc = this->newJc(jcs[jj], resol_);
    jterms_.push_back(jc);
  }
  Log::trace() << "CostFunction: setupTerms Jc added" << std::endl;
  Log::trace() << "CostFunction: setupTerms done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double CostFunction<MODEL>::evaluate(const CtrlVar_ & fguess,
                                     const eckit::Configuration & config,
                                     PostProcessor<State_> post) const {
// Setup terms of cost function
  PostProcessor<State_> pp(post);
  JqTerm_ * jq = jb_->initialize(fguess);
  pp.enrollProcessor(jq);
  for (unsigned jj = 0; jj < jterms_.size(); ++jj) {
    pp.enrollProcessor(jterms_[jj].initialize(fguess));
  }

// Run NL model
  CtrlVar_ mfguess(fguess);
  this->runNL(mfguess, pp);

// Cost function value
  eckit::LocalConfiguration diagnostic;
  if (config.has("diagnostics")) {
    diagnostic = eckit::LocalConfiguration(config, "diagnostics");
  }
  double zzz = jb_->finalize(jq);
  for (unsigned jj = 0; jj < jterms_.size(); ++jj) {
    zzz += jterms_[jj].finalize(diagnostic);
  }
  Log::test() << "CostFunction: Nonlinear J = " << zzz << std::endl;
  return zzz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double CostFunction<MODEL>::linearize(const CtrlVar_ & fguess,
                                      const eckit::Configuration & innerConf,
                                      PostProcessor<State_> post) {
// Read inner loop configuration
  const eckit::LocalConfiguration resConf(innerConf, "resolution");
  const eckit::LocalConfiguration tlmConf(innerConf, "linearmodel");
  eckit::LocalConfiguration diagnostic;
  if (innerConf.has("diagnostics")) {
    diagnostic = eckit::LocalConfiguration(innerConf, "diagnostics");
  }
  const Geometry_ lowres(resConf);

// Setup terms of cost function
  PostProcessor<State_> pp(post);
  JqTerm_ * jq = jb_->initializeTraj(fguess, lowres);
  pp.enrollProcessor(jq);
  for (unsigned jj = 0; jj < jterms_.size(); ++jj) {
    pp.enrollProcessor(jterms_[jj].initializeTraj(fguess, lowres, innerConf));
  }

// Setup linear model (and trajectory)
  tlm_.clear();
  pp.enrollProcessor(new TrajectorySaver<MODEL>(fguess.state()[0], tlmConf, lowres, fguess.modVar(), tlm_));

// Run NL model
  CtrlVar_ mfguess(fguess);
  this->runNL(mfguess, pp);

// Cost function value
  costJb_ = jb_->finalizeTraj(jq);
  costJoJc_ = 0.0;
  for (unsigned jj = 0; jj < jterms_.size(); ++jj) {
    costJoJc_ += jterms_[jj].finalizeTraj(diagnostic);
  }
  double costJ = costJb_ + costJoJc_;
  Log::test() << "CostFunction: Nonlinear J = " << costJ << std::endl;
  return costJ;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostFunction<MODEL>::computeGradientFG(CtrlInc_ & grad) const {
  PostProcessor<Increment_> pp;
  PostProcessorAD<Increment_> costad;
  this->zeroAD(grad);

  for (unsigned jj = 0; jj < jterms_.size(); ++jj) {
    boost::shared_ptr<const GeneralizedDepartures> tmp(jterms_[jj].newGradientFG());
    costad.enrollProcessor(jterms_[jj].setupAD(tmp, grad));
  }

  this->runADJ(grad, costad, pp);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostFunction<MODEL>::addIncrement(CtrlVar_ & xx, const CtrlInc_ & dx,
                                       PostProcessor<Increment_> post) const {
  Log::info() << std::endl;
  Log::info() << "CostFunction::addIncrement: First guess:" << xx << std::endl;
  Log::info() << "CostFunction::addIncrement: Increment:" << dx << std::endl;

  xx.obsVar() += dx.obsVar();
  xx.modVar() += dx.modVar();
  this->addIncr(xx, dx, post);

  Log::info() << "CostFunction::addIncrement: Analysis:" << xx << std::endl;
  Log::test() << "CostFunction::addIncrement: Analysis norm: " << xx.norm() << std::endl;
  Log::info() << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostFunction<MODEL>::resetLinearization() {
  for (unsigned jj = 0; jj < jterms_.size(); ++jj) {
    jterms_[jj].resetLinearization();
  }
  tlm_.clear();
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTFUNCTION_H_
