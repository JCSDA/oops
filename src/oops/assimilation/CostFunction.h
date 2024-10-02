/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2021-2023 UCAR
 * (C) Crown Copyright 2024, the Met Office.
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
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {
  template <typename MODEL, typename OBS> class CostFactory;
  template <typename MODEL, typename OBS> class CostJbState;

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
  typedef CostJo<MODEL, OBS>            CostJo_;
  typedef CostTermBase<MODEL, OBS>      CostBase_;
  typedef Geometry<MODEL>               Geometry_;
  typedef State<MODEL>                  State_;
  typedef Increment<MODEL>              Increment_;

 public:
  CostFunction() = default;
  virtual ~CostFunction() = default;

  virtual double evaluate(CtrlVar_ &, const eckit::Configuration &, PostProcessor<State_>);

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
  virtual void resetLinearization();

/// Compute cost function gradient at first guess (without Jb).
  void computeGradientFG(CtrlInc_ &) const;

/// Append function to call appendObs from jo
  void appendObs(const eckit::Configuration & appendConfig);

/// Access \f$ J_b\f$
  const JbTotal_ & jb() const {return *jb_;}
/// Access terms of the cost function other than \f$ J_b\f$
  const CostJo_ & jo() const {return *jo_;}
  virtual const CostBase_ & jterm(const size_t ii) const {return *(jterms_.at(ii));}
  virtual size_t nterms() const {return jterms_.size();}
  virtual double getCostJoJc() const {return costJoJc_;}

 protected:
  void setupTerms(const eckit::Configuration &);
  std::shared_ptr<JbTotal_> & getNonConstJb() {return jb_;}
  std::shared_ptr<CostJo_> & getNonConstJo() {return jo_;}
  std::vector<std::shared_ptr<CostBase_>> & getJTerms() {return jterms_;}
  virtual const Geometry_ & geometry() const = 0;

 private:
  virtual void addIncr(CtrlVar_ &, const CtrlInc_ &, PostProcessor<Increment_>&) const = 0;

  virtual CostJbState<MODEL, OBS>  * newJb(const eckit::Configuration &,
                                           const Geometry_ &) const = 0;
  virtual CostJo<MODEL, OBS>       * newJo(const eckit::Configuration &) const = 0;
  virtual CostTermBase<MODEL, OBS> * newJc(const eckit::Configuration &,
                                           const Geometry_ &) const = 0;
  virtual void doLinearize(const Geometry_ &, const eckit::Configuration &, CtrlVar_ &, CtrlVar_ &,
                           PostProcessor<State_> &, PostProcessorTLAD<MODEL> &) = 0;
  virtual void finishLinearize() {}

// Data members
  std::shared_ptr<JbTotal_> jb_;
  std::vector<std::shared_ptr<CostBase_>> jterms_;
  std::shared_ptr<CostJo_> jo_;
  std::vector<std::unique_ptr<const Geometry_>> lowres_;

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
    throw std::runtime_error(name + " already registered in cost function factory.");
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
    throw std::runtime_error(id + " does not exist in cost function factory.");
  }
  Log::trace() << "CostFactory::create found cost function type" << std::endl;
  return (*j).second->make(config, comm);
}

// -----------------------------------------------------------------------------
//  Cost Function
// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFunction<MODEL, OBS>::setupTerms(const eckit::Configuration & config) {
  Log::trace() << "CostFunction::setupTerms start" << std::endl;

// Jo
  eckit::LocalConfiguration obsconf(config, "observations");
  jo_.reset(this->newJo(obsconf));
  jterms_.push_back(jo_);
  Log::trace() << "CostFunction::setupTerms Jo added" << std::endl;

// Jb
  CostJbState<MODEL, OBS> * jbs = this->newJb(config, this->geometry());  // constructs background
  jb_.reset(new JbTotal_(jbs, config, this->geometry(), jo_->obspaces()));
  Log::trace() << "CostFunction::setupTerms Jb added" << std::endl;

// Other constraints
  std::vector<eckit::LocalConfiguration> jcs;
  config.get("constraints", jcs);
  for (size_t jj = 0; jj < jcs.size(); ++jj) {
    std::shared_ptr<CostTermBase<MODEL, OBS>> jc(this->newJc(jcs[jj], this->geometry()));
    jterms_.push_back(jc);
  }
  Log::trace() << "CostFunction::setupTerms Jc added" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
double CostFunction<MODEL, OBS>::evaluate(CtrlVar_ & fguess,
                                          const eckit::Configuration & innerConf,
                                          PostProcessor<State_> post) {
  const bool linearize = innerConf.getBool("linearize", false);
  Log::trace() << "CostFunction::evaluate start, linearize = " << linearize << std::endl;

// Setup terms of cost function
  jb_->setPostProc(fguess, innerConf, post);
  for (size_t jj = 0; jj < jterms_.size(); ++jj) {
    jterms_[jj]->setPostProc(fguess, innerConf, post);
  }

  if (linearize) {
//  Inner loop resolution
    const eckit::LocalConfiguration resConf(innerConf, "geometry");
    lowres_.push_back(std::make_unique<Geometry_>(resConf, this->geometry().getComm(),
                                                  this->geometry().timeComm()));

//  Setup trajectory for terms of cost function
    PostProcessorTLAD<MODEL> pptraj;
    jb_->setPostProcTraj(fguess, innerConf, *lowres_.back(), pptraj);
    for (size_t jj = 0; jj < jterms_.size(); ++jj) {
      jterms_[jj]->setPostProcTraj(fguess, innerConf, *lowres_.back(), pptraj);
    }

//  Setup specific linearization if needed (including TLM)
    this->doLinearize(*lowres_.back(), innerConf, jb_->getBackground(), fguess, post, pptraj);
  }

// Run NL model
  CtrlVar_ xxtmp(fguess);
  this->runNL(xxtmp, post);

// Cost function value (except Jb)
  costJoJc_ = 0.0;
  for (size_t jj = 0; jj < jterms_.size(); ++jj) {
    costJoJc_ += jterms_[jj]->computeCost();  // needed before jb_->computeCostTraj()
  }

  if (linearize) {
//  Finalize trajectory setup
    this->finishLinearize();  // Used only for FGAT
    jb_->computeCostTraj();   // constructs B (requires H(x)/QC for VarBC part of B)
    for (size_t jj = 0; jj < jterms_.size(); ++jj) {
      jterms_[jj]->computeCostTraj();
    }
  }

// Cost function value (Jb)
  costJb_ = jb_->computeCost();  // Should be after linearization because requires B

// Print cost function for test (hack so we can update references separately)
  for (size_t jj = 0; jj < jterms_.size(); ++jj) {
    jterms_[jj]->printCostTestHack();
  }

  double zzz = costJb_ + costJoJc_;
  Log::info() << "CostFunction: Nonlinear J = " << zzz << std::endl;
  Log::test() << "CostFunction: Nonlinear J = " << zzz << std::endl;

//  First-guess gradient of Jb needs to be reset to zero before minimisation
//  when the Pert members of Control-Pert EDA are run
  if (innerConf.has("control pert")) this->getNonConstJb()->zeroGradientFG();

  Log::trace() << "CostFunction::evaluate done" << std::endl;
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
    std::shared_ptr<const GeneralizedDepartures> tmp(jterms_[jj]->newGradientFG());
    jterms_[jj]->computeCostAD(tmp, grad, costad);
  }

  this->runADJ(grad, costad, pp);

  for (size_t jj = 0; jj < jterms_.size(); ++jj) {
    jterms_[jj]->setPostProcAD();
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
    jterms_[jj]->resetLinearization();
  }
  Log::trace() << "CostFunction::resetLinearization done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFunction<MODEL, OBS>::appendObs(const eckit::Configuration & appendConfig) {
  Log::trace() << "CostFunction::appendObs start" << std::endl;
  jo_->appendObs(appendConfig);
  Log::trace() << "CostFunction::appendObs done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTFUNCTION_H_
