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

#include <boost/noncopyable.hpp>
#include <map>
#include <string>

#include "eckit/config/Configuration.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualVector.h"
#include "oops/assimilation/HMatrix.h"
#include "oops/assimilation/HtMatrix.h"
#include "oops/assimilation/Minimizer.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"

#include "util/abor1_cpp.h"
#include "util/dot_product.h"
#include "util/formats.h"
#include "util/Logger.h"
#include "util/PrintAdjTest.h"

namespace oops {

// -----------------------------------------------------------------------------

/// A Minimizer knows how to minimize a cost function
template<typename MODEL> class Minimizer : private boost::noncopyable {
  typedef CostFunction<MODEL>      CostFct_;
  typedef ControlIncrement<MODEL>  CtrlInc_;
  typedef DualVector<MODEL>        Dual_;
  typedef HMatrix<MODEL>           H_;
  typedef HtMatrix<MODEL>          Ht_;
  typedef State<MODEL>             State_;

 public:
  explicit Minimizer(const CostFct_ & J): J_(J), outerIteration_(0) {};
  virtual ~Minimizer() {}
  ControlIncrement<MODEL> * minimize(const eckit::Configuration &);
  virtual const std::string classname() const =0;

 private:
  virtual ControlIncrement<MODEL> * doMinimize(const eckit::Configuration &) =0;

  void adjTests(const eckit::Configuration &);
  void adjModelTest(const Ht_ &, const H_ &);
  void adjObsTest(const Ht_ &, const H_ &);

  void tlmTests(const eckit::Configuration &);
  void tlmApproxTest(const H_ &);
  void tlmTaylorTest(const H_ &);

  void tlmPropagTest(const eckit::Configuration & config, const CtrlInc_ &);

  const CostFct_ & J_;
  int outerIteration_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
ControlIncrement<MODEL> * Minimizer<MODEL>::minimize(const eckit::Configuration & config) {
  // TLM tests
  this->tlmTests(config);  

  // ADJ tests
  this->adjTests(config);

  // Minimize
  ControlIncrement<MODEL> * dx = this->doMinimize(config);

  // TLM propagation test
  this->tlmPropagTest(config, *dx);  

  // Update outer loop counter
  outerIteration_++; 

  return dx;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Minimizer<MODEL>::tlmTests(const eckit::Configuration & config) {
// Tangent Linear tests

  if (config.has("onlineDiagnostics")) {
    const eckit::LocalConfiguration onlineDiag(config, "onlineDiagnostics");
    bool runTlmTaylorTest = onlineDiag.getBool("tlmTaylorTest");
    bool runTlmApproxTest = onlineDiag.getBool("tlmApproxTest");

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

template<typename MODEL>
void Minimizer<MODEL>::adjTests(const eckit::Configuration & config) {
// Adjoint tests

  if (config.has("onlineDiagnostics")) {
    const eckit::LocalConfiguration onlineDiag(config, "onlineDiagnostics");
    bool runAdjTlmTest    = onlineDiag.getBool("adjTlmTest");
    bool runAdjObsTest    = onlineDiag.getBool("adjObsTest");

    const H_  H(J_);
    const Ht_ Ht(J_);
  
    // Online tests
    // Model adjoint test
    if (runAdjTlmTest) this->adjModelTest(Ht, H);
  
    // Obs adjoint test
    if (runAdjObsTest) this->adjObsTest(Ht, H);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Minimizer<MODEL>::tlmApproxTest(const H_ & H) {
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
  ControlVariable<MODEL> mxx(J_.jb().getBackground());
  J_.runNL(mxx, pp);

// run NL for perturbed initial condition
  ControlVariable<MODEL> mpertxx(J_.jb().getBackground());
  J_.addIncrement(mpertxx, dx);
  J_.runNL(mpertxx, pp);

// run TL for the perturbation
  Dual_ dummy;
  CtrlInc_ mdx(dx);
  H.multiply(mdx, dummy);

// print log
  Log::info() << std::endl
              << "TLM Linear Approximation Test: done." 
              << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Minimizer<MODEL>::tlmPropagTest(const eckit::Configuration & config, 
                                     const CtrlInc_ & dx) {
/* TL propagation test:
   calculate and store:
     M'(dx)
   where dx comes from minimization stopped at requested iteration 
   
   TO DO:
     write to file depending on requriements */

  if (config.has("onlineDiagnostics")) {
    const eckit::LocalConfiguration onlineDiag(config, "onlineDiagnostics");
    bool runTlmPropagTest = onlineDiag.getBool("tlmPropagTest");

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
      Log::info() << std::endl
                  << "TLM Propagation Test: done." 
                  << std::endl;
    }
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Minimizer<MODEL>::tlmTaylorTest(const H_ & H) {
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
  ControlVariable<MODEL> mxx(J_.jb().getBackground());
  J_.runNL(mxx, pp);

// run TL for the un-damped perturbation: M'(dx)
  Dual_ dummy;
  CtrlInc_ mdx(dx);
  H.multiply(mdx, dummy);

// loop over decreasing increments
  for (unsigned int i = 0; i < 14; ++i) {
     // ||p M'(dx)||
     CtrlInc_ pmdx(mdx);
     pmdx *= 1./pow(10.0,i);
     double denom = sqrt(dot_product(pmdx, pmdx));

     // run perturbed NL: M(x+pdx)
     ControlVariable<MODEL> mpertxx(J_.jb().getBackground());
     CtrlInc_ pdx(dx);
     pdx *= 1./pow(10.0,i);
     J_.addIncrement(mpertxx, pdx);
     J_.runNL(mpertxx, pp);

     // ||M(x+pdx) - M(x)||
     CtrlInc_ diff_nl(mdx, false);
     diff_nl.state()[0].diff(mpertxx.state()[0], mxx.state()[0]);
     double nom = sqrt(dot_product(diff_nl, diff_nl));

     // print results
     Log::info() << std::endl
                 << "TLM Taylor test:  p = " << std::setw(8) << 1./pow(10.0,i)
                 << ", ||M(x) - M(x+p dx)|| / ||p M'(dx)|| = "
                 << util::full_precision(1.+std::abs(1.-nom/denom))
                 << std::endl << std::endl;
  }
  Log::info() << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Minimizer<MODEL>::adjModelTest(const Ht_ & Ht,
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

// run TL
  Dual_ dummy;
  CtrlInc_ mdx1(dx1);
  H.multiply(mdx1, dummy, false);

// run ADJ
  dummy.zero(); 
  CtrlInc_ mtdx2(dx2);
  mtdx2.state()[0].updateTime(mdx1.state()[0].validTime() - dx1.state()[0].validTime());
  Ht.multiply(dummy, mtdx2, false);

// calculate FWD < M dx1, dx2 >
  double adj_tst_fwd = dot_product(mdx1, dx2);

// calculate BWD < dx1, Mt dx2 >
  double adj_tst_bwd = dot_product(dx1, mtdx2);

// print results
  Log::info() << "Model Adjoint Test: " << outerIteration_ << std::endl
              << util::PrintAdjTest(adj_tst_fwd, adj_tst_bwd, "M")
              << std::endl << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Minimizer<MODEL>::adjObsTest(const Ht_ & Ht,
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
  J_.jb().randomize(dx2);

// run TL
  H.multiply(dx1, hdx1, true);
  H.multiply(dx2, hdx2, true);

// run ADJ
  J_.zeroAD(hthdx2);
  Ht.multiply(hdx2, hthdx2, true);

// calculate FWD < H dx1, hdx2 >
  double adj_tst_fwd = dot_product(hdx1, hdx2);

// calculate BWD < dx1, Ht hdx2 >
  double adj_tst_bwd = dot_product(dx1, hthdx2);

// print results
  Log::info() << "Obs Adjoint Test: " << outerIteration_ << std::endl
              << util::PrintAdjTest(adj_tst_fwd, adj_tst_bwd, "H")
              << std::endl << std::endl;
}

// -----------------------------------------------------------------------------

/// Minimizer Factory
template <typename MODEL>
class MinFactory {
  typedef CostFunction<MODEL>        CostFct_;
 public:
  static Minimizer<MODEL> * create(const eckit::Configuration &, const CostFct_ &);
  virtual ~MinFactory() { getMakers().clear(); }
 protected:
  explicit MinFactory(const std::string &);
 private:
  virtual Minimizer<MODEL> * make(const eckit::Configuration &, const CostFct_ &) =0;
  static std::map < std::string, MinFactory<MODEL> * > & getMakers() {
    static std::map < std::string, MinFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class FCT>
class MinMaker : public MinFactory<MODEL> {
  typedef CostFunction<MODEL>        CostFct_;
  virtual Minimizer<MODEL> * make(const eckit::Configuration & conf, const CostFct_ & J) {
    return new FCT(conf, J);
  }
 public:
  explicit MinMaker(const std::string & name) : MinFactory<MODEL>(name) {}
};

// =============================================================================

template <typename MODEL>
MinFactory<MODEL>::MinFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    Log::error() << name << " already registered in minimizer factory." << std::endl;
    ABORT("Element already registered in MinFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Minimizer<MODEL>* MinFactory<MODEL>::create(const eckit::Configuration & config,
                                            const CostFct_ & J) {
  std::string id = config.getString("algorithm");
  Log::info() << "Minimizer algorithm=" << id << std::endl;
  typename std::map<std::string, MinFactory<MODEL>*>::iterator j = getMakers().find(id);
  if (j == getMakers().end()) {
    Log::error() << id << " does not exist in minimizer factory." << std::endl;
    ABORT("Element does not exist in MinFactory.");
  }
  return (*j).second->make(config, J);
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_MINIMIZER_H_
