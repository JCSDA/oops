/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_GENERIC_ERRORCOVARIANCEBUMP_H_
#define OOPS_GENERIC_ERRORCOVARIANCEBUMP_H_

#include <string>
#include <vector>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/assimilation/GMRESR.h"
#include "oops/assimilation/Increment4D.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
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

/// Model space error covariance on generic unstructured grid

template <typename MODEL>
class ErrorCovarianceBUMP : public oops::ModelSpaceCovarianceBase<MODEL>,
                            public util::Printable,
                            private util::ObjectCounter<ErrorCovarianceBUMP<MODEL>>,
                            private boost::noncopyable {
  typedef Geometry<MODEL>                         Geometry_;
  typedef Increment<MODEL>                        Increment_;
  typedef Increment4D<MODEL>                      Increment4D_;
  typedef State<MODEL>                            State_;
  typedef ParametersBUMP<MODEL>                   Parameters_;
  typedef StateEnsemble<MODEL>                    Ensemble_;
  typedef boost::shared_ptr<StateEnsemble<MODEL>> EnsemblePtr_;

 public:
  static const std::string classname() {return "oops::ErrorCovarianceBUMP";}

  ErrorCovarianceBUMP(const Geometry_ &, const Variables &,
                      const eckit::Configuration &, const State_ &, const State_ &);
  virtual ~ErrorCovarianceBUMP();

  void randomize(Increment_ &) const;

 private:
  void doRandomize(Increment_ &) const override;
  void doMultiply(const Increment_ &, Increment_ &) const override;
  void doInverseMultiply(const Increment_ &, Increment_ &) const override;

  void print(std::ostream &) const override;

  int keyBUMP_;
};

// =============================================================================

template<typename MODEL>
ErrorCovarianceBUMP<MODEL>::ErrorCovarianceBUMP(const Geometry_ & resol,
                                                const Variables & vars,
                                                const eckit::Configuration & conf,
                                                const State_ & xb, const State_ & fg)
  : ModelSpaceCovarianceBase<MODEL>(xb, fg, resol, conf), keyBUMP_(0)
{
  Log::trace() << "ErrorCovarianceBUMP::ErrorCovarianceBUMP starting" << std::endl;

// Setup timeslots
  std::vector<util::DateTime> timeslots;
  timeslots.push_back(xb.validTime());

// Setup ensemble of perturbations
  EnsemblePtr_ ens(new Ensemble_());

// Setup pseudo ensemble
  EnsemblePtr_ pseudo_ens(new Ensemble_());

// Setup parameters
  Parameters_ param(resol, vars, timeslots, ens, pseudo_ens, conf);

// Get key
  keyBUMP_ = param.get_bump();

  Log::trace() << "ErrorCovarianceBUMP::ErrorCovarianceBUMP done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovarianceBUMP<MODEL>::~ErrorCovarianceBUMP() {
  Log::trace() << "ErrorCovarianceBUMP<MODEL>::~ErrorCovarianceBUMP starting" << std::endl;
  util::Timer timer(classname(), "~ErrorCovarianceBUMP");
  delete_oobump_f90(keyBUMP_);
  Log::trace() << "ErrorCovarianceBUMP<MODEL>::~ErrorCovarianceBUMP done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovarianceBUMP<MODEL>::randomize(Increment_ & dx) const {
  Log::trace() << "ErrorCovarianceBUMP<MODEL>::randomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");
  int colocated;
  get_oobump_colocated_f90(keyBUMP_, colocated);
  UnstructuredGrid ug(colocated);
  dx.ug_coord(ug);
  randomize_oobump_nicas_f90(keyBUMP_, ug.toFortran());
  dx.field_from_ug(ug);
  Log::trace() << "ErrorCovarianceBUMP<MODEL>::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovarianceBUMP<MODEL>::doMultiply(const Increment_ & dxi,
                                            Increment_ & dxo) const {
  Log::trace() << "ErrorCovarianceBUMP<MODEL>::doMultiply starting" << std::endl;
  util::Timer timer(classname(), "doMultiply");
  int colocated;
  get_oobump_colocated_f90(keyBUMP_, colocated);
  UnstructuredGrid ug(colocated);
  dxi.field_to_ug(ug);
  multiply_oobump_nicas_f90(keyBUMP_, ug.toFortran());
  dxo.field_from_ug(ug);
  Log::trace() << "ErrorCovarianceBUMP<MODEL>::doMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovarianceBUMP<MODEL>::doInverseMultiply(const Increment_ & dxi,
                                                   Increment_ & dxo) const {
  Log::trace() << "ErrorCovarianceBUMP<MODEL>::doInverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "doInverseMultiply");
  IdentityMatrix<Increment_> Id;
  dxo.zero();
  GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);
  Log::trace() << "ErrorCovarianceBUMP<MODEL>::doInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovarianceBUMP<MODEL>::doRandomize(Increment_ & dx) const {
  Log::trace() << "ErrorCovarianceBUMP<MODEL>::doRandomize starting" << std::endl;
  util::Timer timer(classname(), "doRandomize");
  dx.random();
  Log::trace() << "ErrorCovarianceBUMP<MODEL>::doRandomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovarianceBUMP<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ErrorCovarianceBUMP<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << "ErrorCovarianceBUMP<MODEL>::print not implemented";
  Log::trace() << "ErrorCovarianceBUMP<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_ERRORCOVARIANCEBUMP_H_
