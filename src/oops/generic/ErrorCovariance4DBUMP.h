/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_GENERIC_ERRORCOVARIANCE4DBUMP_H_
#define OOPS_GENERIC_ERRORCOVARIANCE4DBUMP_H_

#include <string>
#include <vector>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/assimilation/GMRESR.h"
#include "oops/assimilation/Increment4D.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/ModelSpaceCovariance4DBase.h"
#include "oops/base/Variables.h"
#include "oops/generic/oobump_f.h"
#include "oops/generic/ParametersBUMP.h"
#include "oops/generic/UnstructuredGrid.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class LocalConfiguration;
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

/// Model space error 4D covariance on generic unstructured grid

template <typename MODEL>
class ErrorCovariance4DBUMP : public oops::ModelSpaceCovariance4DBase<MODEL>,
                              public util::Printable,
                              private util::ObjectCounter<ErrorCovariance4DBUMP<MODEL> >,
                              private boost::noncopyable {
  typedef Geometry<MODEL>                         Geometry_;
  typedef Increment4D<MODEL>                      Increment4D_;
  typedef State<MODEL>                            State_;
  typedef State4D<MODEL>                          State4D_;
  typedef ParametersBUMP<MODEL>                   Parameters_;
  typedef StateEnsemble<MODEL>                    Ensemble_;
  typedef boost::shared_ptr<StateEnsemble<MODEL>> EnsemblePtr_;

 public:
  static const std::string classname() {return "oops::ErrorCovariance4DBUMP";}

  ErrorCovariance4DBUMP(const Geometry_ &, const Variables &,
                        const eckit::Configuration &, const State4D_ &, const State4D_ &);
  virtual ~ErrorCovariance4DBUMP();

 private:
  void doRandomize(Increment4D_ &) const override;
  void doMultiply(const Increment4D_ &, Increment4D_ &) const override;
  void doInverseMultiply(const Increment4D_ &, Increment4D_ &) const override;

  void print(std::ostream &) const override;
  std::vector<util::DateTime> timeslots_;
  int keyBUMP_;
};

// =============================================================================

template<typename MODEL>
ErrorCovariance4DBUMP<MODEL>::ErrorCovariance4DBUMP(const Geometry_ & resol,
                                                    const Variables & vars,
                                                    const eckit::Configuration & conf,
                                                    const State4D_ & xb, const State4D_ & fg)
  : ModelSpaceCovariance4DBase<MODEL>(xb, fg, resol, conf), keyBUMP_(0)
{
  Log::trace() << "ErrorCovariance4DBUMP::ErrorCovariance4DBUMP starting" << std::endl;

//  Setup timeslots
  for (unsigned jsub = 0; jsub < xb.size(); ++jsub) {
    timeslots_.push_back(xb[jsub].validTime());
  }

// Setup ensemble of perturbations
  EnsemblePtr_ ens(new Ensemble_());

// Setup pseudo ensemble
  EnsemblePtr_ pseudo_ens(new Ensemble_());

// Setup parameters
  Parameters_ param(resol, vars, timeslots_, ens, pseudo_ens, conf);

// Get key
  keyBUMP_ = param.get_bump();

  Log::trace() << "ErrorCovariance4DBUMP::ErrorCovariance4DBUMP done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance4DBUMP<MODEL>::~ErrorCovariance4DBUMP() {
  Log::trace() << "ErrorCovariance4DBUMP<MODEL>::~ErrorCovariance4DBUMP starting" << std::endl;
  util::Timer timer(classname(), "~ErrorCovariance4DBUMP");
  delete_oobump_f90(keyBUMP_);
  Log::trace() << "ErrorCovariance4DBUMP<MODEL>::~ErrorCovariance4DBUMP done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance4DBUMP<MODEL>::doRandomize(Increment4D_ & dx) const {
  Log::trace() << "ErrorCovariance4DBUMP<MODEL>::doRandomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");
  int colocated;
  get_oobump_colocated_f90(keyBUMP_, colocated);
  UnstructuredGrid ug(colocated, timeslots_.size());
  dx.ug_coord(ug);
  randomize_oobump_nicas_f90(keyBUMP_, ug.toFortran());
  dx.field_from_ug(ug);
  Log::trace() << "ErrorCovariance4DBUMP<MODEL>::doRandomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance4DBUMP<MODEL>::doMultiply(const Increment4D_ & dxi,
                                              Increment4D_ & dxo) const {
  Log::trace() << "ErrorCovariance4DBUMP<MODEL>::doMultiply starting" << std::endl;
  util::Timer timer(classname(), "doMultiply");
  int colocated;
  get_oobump_colocated_f90(keyBUMP_, colocated);
  UnstructuredGrid ug(colocated, timeslots_.size());
  dxi.field_to_ug(ug);
  multiply_oobump_nicas_f90(keyBUMP_, ug.toFortran());
  dxo.field_from_ug(ug);
  Log::trace() << "ErrorCovariance4DBUMP<MODEL>::doMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance4DBUMP<MODEL>::doInverseMultiply(const Increment4D_ & dxi,
                                                     Increment4D_ & dxo) const {
  Log::trace() << "ErrorCovariance4DBUMP<MODEL>::doInverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "doInverseMultiply");
  IdentityMatrix<Increment4D_> Id;
  dxo.zero();
  GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);
  Log::trace() << "ErrorCovariance4DBUMP<MODEL>::doInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance4DBUMP<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ErrorCovariance4DBUMP<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << "ErrorCovariance4DBUMP<MODEL>::print not implemented";
  Log::trace() << "ErrorCovariance4DBUMP<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_ERRORCOVARIANCE4DBUMP_H_