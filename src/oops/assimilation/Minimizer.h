/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_MINIMIZER_H_
#define OOPS_ASSIMILATION_MINIMIZER_H_

#include <map>
#include <string>
#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualVector.h"
#include "oops/assimilation/HMatrix.h"
#include "oops/assimilation/HtMatrix.h"
#include "oops/assimilation/MinimizerUtils.h"
#include "oops/base/State.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"
#include "oops/util/PrintAdjTest.h"

namespace oops {

// -----------------------------------------------------------------------------

/// A Minimizer knows how to minimize a cost function
template<typename MODEL, typename OBS> class Minimizer : private boost::noncopyable {
  typedef CostFunction<MODEL, OBS>      CostFct_;
  typedef ControlIncrement<MODEL, OBS>  CtrlInc_;
  typedef DualVector<MODEL, OBS>        Dual_;
  typedef HMatrix<MODEL, OBS>           H_;
  typedef HtMatrix<MODEL, OBS>          Ht_;
  typedef State<MODEL>                  State_;

 public:
  explicit Minimizer(const CostFct_ & J): J_(J), outerIteration_(0) {}
  virtual ~Minimizer() {}
  ControlIncrement<MODEL, OBS> * minimize(const eckit::Configuration &);
  virtual const std::string classname() const = 0;

 private:
  virtual ControlIncrement<MODEL, OBS> * doMinimize(const eckit::Configuration &) = 0;

  void adjTests(const eckit::Configuration &);
  void adjModelTest(const Ht_ &, const H_ &);
  void adjObsTest(const Ht_ &, const H_ &);
  void adjGeneralizedObsOpTest(const Ht_ &, const H_ &);

  void tlmTests(const eckit::Configuration &);
  void tlmApproxTest(const H_ &);
  void tlmTaylorTest(const H_ &);

  void tlmPropagTest(const eckit::Configuration & config, const CtrlInc_ &);
  void outputAMinusB(const eckit::Configuration & config, CtrlInc_ controlInc);

  const CostFct_ & J_;
  int outerIteration_;
};

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
ControlIncrement<MODEL, OBS> *
Minimizer<MODEL, OBS>::minimize(const eckit::Configuration & config) {
  // TLM tests
  this->tlmTests(config);

  // ADJ tests
  this->adjTests(config);

  // Minimize
  ControlIncrement<MODEL, OBS> * dx = this->doMinimize(config);
  this->outputAMinusB(config, *dx);

  // Write increment
  writeIncrement(config, *dx, outerIteration_);

  // TLM propagation test
  this->tlmPropagTest(config, *dx);

  // Update outer loop counter
  outerIteration_++;

  return dx;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void Minimizer<MODEL, OBS>::tlmTests(const eckit::Configuration & config) {
// Tangent Linear tests

  if (config.has("online diagnostics")) {
    const eckit::LocalConfiguration onlineDiag(config, "online diagnostics");
    bool runTlmTaylorTest = onlineDiag.getBool("tlm taylor test", false);
    bool runTlmApproxTest = onlineDiag.getBool("tlm approx test", false);

    const H_  H(J_);
    const Ht_ Ht(J_);

    // Online tests
    // TLM linear approximation test
    if (outerIteration_ == 0 && runTlmApproxTest) this->tlmApproxTest(H);

    // TLM Taylor test
    if (outerIteration_ == 0 && runTlmTaylorTest) this->tlmTaylorTest(H);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void Minimizer<MODEL, OBS>::adjTests(const eckit::Configuration & config) {
// Adjoint tests

  if (config.has("online diagnostics")) {
    const eckit::LocalConfiguration onlineDiag(config, "online diagnostics");
    bool runAdjTlmTest    = onlineDiag.getBool("adj tlm test", false);
    bool runAdjObsTest    = onlineDiag.getBool("adj obs test", false);
    bool runAdjGeneralizedObsOpTest = onlineDiag.getBool("adj generalized obs op test", false);

    const H_  H(J_);
    const Ht_ Ht(J_);

    // Online tests
    // Model adjoint test
    if (runAdjTlmTest) this->adjModelTest(Ht, H);

    // Obs adjoint test
    if (runAdjObsTest) this->adjObsTest(Ht, H);

    if (runAdjGeneralizedObsOpTest) this->adjGeneralizedObsOpTest(Ht, H);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void Minimizer<MODEL, OBS>::tlmApproxTest(const H_ & H) {
/* TL approx test:
   calculate and store:
     M(x + dx), M(x), M'(dx)
   where dx is initialized using randomization

   TO DO:
     write to file depending on requriements */

  Log::info() << std::endl
              << "TLM Linear Approximation Test - starting: " << outerIteration_
              << std::endl << std::endl;

// initialize incremets
  CtrlInc_ dx(J_.jb());

// randomize increments
  J_.jb().randomize(dx);

// run NL for unperturbed initial condition
  PostProcessor<State_> pp;
  ControlVariable<MODEL, OBS> mxx(J_.jb().getBackground());
  J_.runNL(mxx, pp);

// run NL for perturbed initial condition
  ControlVariable<MODEL, OBS> mpertxx(J_.jb().getBackground());
  J_.addIncrement(mpertxx, dx);
  J_.runNL(mpertxx, pp);

// run TL for the perturbation
  Dual_ dummy;
  CtrlInc_ mdx(dx);
  H.multiply(mdx, dummy);

// print log
  Log::info() << std::endl << "TLM Linear Approximation Test: done." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void Minimizer<MODEL, OBS>::tlmPropagTest(const eckit::Configuration & config,
                                     const CtrlInc_ & dx) {
/* TL propagation test:
   calculate and store:
     M'(dx)
   where dx comes from minimization stopped at requested iteration

   TO DO:
     write to file depending on requriements */

  if (config.has("online diagnostics")) {
    const eckit::LocalConfiguration onlineDiag(config, "online diagnostics");
    bool runTlmPropagTest = onlineDiag.getBool("tlm propagation test", false);

    if (runTlmPropagTest) {
      // print log
      Log::info() << "TLM Propagation Test - starting: " << outerIteration_
                  << std::endl << std::endl;

      // construct propagation matrix
      const H_ H(J_);

      // run TL for the input perturbation: M'(dx)
      Dual_ dummy;
      CtrlInc_ mdx(dx);
      H.multiply(mdx, dummy);

      // print log
      Log::info() << std::endl << "TLM Propagation Test: done." << std::endl;
    }
  }
}

template<typename MODEL, typename OBS>
void Minimizer<MODEL, OBS>::outputAMinusB(const eckit::Configuration & config,
                                          CtrlInc_ controlInc) {
  if (config.getBool("output a minus b", false)) {
    const auto H = HMatrix<MODEL, OBS>(J_);
    DualVector<MODEL, OBS> aMinusB;
    H.multiply(controlInc, aMinusB);
    const std::string depName = "output_a_minus_b" + std::to_string(outerIteration_);
    aMinusB.saveDep(depName);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void Minimizer<MODEL, OBS>::tlmTaylorTest(const H_ & H) {
/* TL Taylor test:
   ||M(x + p dx) - M(x)|| / ||M'(p dx)|| -> 1
   where p is a scalar damping factor and dx is initialized
   using randomization */

  Log::info() << std::endl
              << "TLM Taylor Test: " << outerIteration_
              << std::endl;

// initialize incremets
  CtrlInc_ dx(J_.jb());

// randomize increments
  J_.jb().randomize(dx);

// run NL for unperturbed initial condition: M(x)
  PostProcessor<State_> pp;
  ControlVariable<MODEL, OBS> mxx(J_.jb().getBackground());
  J_.runNL(mxx, pp);

// run TL for the un-damped perturbation: M'(dx)
  Dual_ dummy;
  CtrlInc_ mdx(dx);
  H.multiply(mdx, dummy);

// loop over decreasing increments
  for (unsigned int jj = 0; jj < 14; ++jj) {
     // ||p M'(dx)||
     CtrlInc_ pmdx(mdx);
     pmdx *= 1./pow(10.0, jj);
     double denom = sqrt(dot_product(pmdx, pmdx));

     // run perturbed NL: M(x+pdx)
     ControlVariable<MODEL, OBS> mpertxx(J_.jb().getBackground());
     CtrlInc_ pdx(dx);
     pdx *= 1./pow(10.0, jj);
     J_.addIncrement(mpertxx, pdx);
     J_.runNL(mpertxx, pp);

     // ||M(x+pdx) - M(x)||
     CtrlInc_ diff_nl(mdx, false);
     diff_nl.state().diff(mpertxx.state(), mxx.state());
     double nom = sqrt(dot_product(diff_nl, diff_nl));

     // print results
     Log::info() << std::endl
                 << "TLM Taylor test:  p = " << std::setw(8) << 1./pow(10.0, jj)
                 << ", ||M(x) - M(x+p dx)|| / ||p M'(dx)|| = "
                 << util::full_precision(1.+std::abs(1.-nom/denom))
                 << std::endl;
  }
  Log::info() << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void Minimizer<MODEL, OBS>::adjModelTest(const Ht_ & Ht,
                                    const H_ & H) {
/* Perform the adjoint test of the linear model
   <M dx1, dx2> - <dx1, Mt dx2>)/<M dx1, dx2>
   where dx1 and dx2 are increments obtained through randomization */

// Model adjoint test
  CtrlInc_ dx1(J_.jb());
  CtrlInc_ dx2(J_.jb());

// randomize increments
  J_.jb().randomize(dx1);
  J_.jb().randomize(dx2);

/* zero out modvar and obsvar
   since only model is being tested
*/
  dx1.obsVar().zero();
  dx1.modVar().zero();
  dx2.obsVar().zero();
  dx2.modVar().zero();

// run TL
  Dual_ dummy;
  CtrlInc_ mdx1(dx1);
  H.multiply(mdx1, dummy, false);

// run ADJ
  dummy.zero();
  CtrlInc_ mtdx2(dx2);
  mtdx2.state().updateTime(mdx1.state().validTime() - dx1.state().validTime());
  Ht.multiply(dummy, mtdx2, false);

// calculate FWD < M dx1, dx2 >
  dx2.state().updateTime(mdx1.state().validTime() - dx1.state().validTime());
  double adj_tst_fwd = dot_product(mdx1, dx2);

// calculate BWD < dx1, Mt dx2 >
  double adj_tst_bwd = dot_product(dx1, mtdx2);

// print results
  Log::info() << "Model Adjoint Test: " << outerIteration_ << std::endl
              << util::PrintAdjTest(adj_tst_fwd, adj_tst_bwd, "M")
              << std::endl << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void Minimizer<MODEL, OBS>::adjObsTest(const Ht_ & Ht,
                                  const H_ & H) {
/* Perform the adjoint test of the linear observation operator
   (<H dx, dy> - <dx, Ht dy>)/<H dx, dy>
   where dx is an increment obtained through randomization and
   dy is a DualVector obtained by first creating an increment
   through randomization and then applying a linear observation
   operator */

// Model adjoint test
  CtrlInc_ dx1(J_.jb());
  CtrlInc_ dx2(J_.jb());
  CtrlInc_ hthdx2(J_.jb());
  Dual_ hdx1;
  Dual_ hdx2;

// randomize increments
  J_.jb().randomize(dx1);
  CtrlInc_ dx1Copy(dx1);
  J_.jb().randomize(dx2);

// run TL
  H.multiply(dx1, hdx1, true);
  H.multiply(dx2, hdx2, true);
  Dual_ hdx2Copy(hdx2);

// run ADJ
  J_.zeroAD(hthdx2);
  Ht.multiply(hdx2, hthdx2, true);

// calculate FWD < H dx1, hdx2 >
  double adj_tst_fwd = dot_product(hdx1, hdx2Copy);

// calculate BWD < dx1, Ht hdx2 >
  double adj_tst_bwd = dot_product(dx1Copy, hthdx2);

// print results
  Log::info() << "Obs Adjoint Test: " << outerIteration_ << std::endl
              << util::PrintAdjTest(adj_tst_fwd, adj_tst_bwd, "H")
              << std::endl << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void Minimizer<MODEL, OBS>::adjGeneralizedObsOpTest(const Ht_ & Ht, const H_ & H) {
  // Random vector in model space at initial time
  CtrlInc_ dxCopy(J_.jb());
  J_.jb().randomize(dxCopy);

  // Random vector in observation space at initial time
  CtrlInc_ dxDummy(J_.jb());
  J_.jb().randomize(dxDummy);
  Dual_ dyCopy;
  H.multiply(dxDummy, dyCopy, false);

  // Initialize terms of adjoint test of G
  CtrlInc_ dx(dxCopy);
  Dual_ Gdx;
  Dual_ dy(dyCopy);
  CtrlInc_ Gtdy(J_.jb());
  Gtdy.state().updateTime(dxDummy.state().validTime() - dx.state().validTime());;

  // Compute terms
  H.multiply(dxCopy, Gdx, false);
  Ht.multiply(dyCopy, Gtdy, false);
  double adj_tst_fwd = dot_product(Gdx, dy);
  double adj_tst_bwd = dot_product(dx, Gtdy);

  Log::info() << "Generalized Obs Op Adjoint Test: " << outerIteration_ << std::endl
              << util::PrintAdjTest(adj_tst_fwd, adj_tst_bwd, "G")
              << std::endl << std::endl;
}

// -----------------------------------------------------------------------------

/// Minimizer Factory
template <typename MODEL, typename OBS>
class MinFactory {
  typedef CostFunction<MODEL, OBS>        CostFct_;
 public:
  static Minimizer<MODEL, OBS> * create(const eckit::Configuration &, const CostFct_ &);
  virtual ~MinFactory() = default;
 protected:
  explicit MinFactory(const std::string &);
 private:
  virtual Minimizer<MODEL, OBS> * make(const eckit::Configuration &, const CostFct_ &) = 0;
  static std::map < std::string, MinFactory<MODEL, OBS> * > & getMakers() {
    static std::map < std::string, MinFactory<MODEL, OBS> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class OBS, class FCT>
class MinMaker : public MinFactory<MODEL, OBS> {
  typedef CostFunction<MODEL, OBS>        CostFct_;
  virtual Minimizer<MODEL, OBS> * make(const eckit::Configuration & conf, const CostFct_ & J) {
    return new FCT(conf, J);
  }
 public:
  explicit MinMaker(const std::string & name) : MinFactory<MODEL, OBS>(name) {}
};

// =============================================================================

template <typename MODEL, typename OBS>
MinFactory<MODEL, OBS>::MinFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in minimizer factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
Minimizer<MODEL, OBS>* MinFactory<MODEL, OBS>::create(const eckit::Configuration & config,
                                            const CostFct_ & J) {
  std::string id = config.getString("algorithm");
  Log::info() << "Minimizer algorithm=" << id << std::endl;
  typename std::map<std::string, MinFactory<MODEL, OBS>*>::iterator j = getMakers().find(id);
  if (j == getMakers().end()) {
    throw std::runtime_error(id + " does not exist in minimizer factory.");
  }
  return (*j).second->make(config, J);
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_MINIMIZER_H_
