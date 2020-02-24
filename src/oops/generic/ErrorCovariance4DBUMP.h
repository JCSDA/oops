/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_ERRORCOVARIANCE4DBUMP_H_
#define OOPS_GENERIC_ERRORCOVARIANCE4DBUMP_H_

#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/assimilation/GMRESR.h"
#include "oops/assimilation/Increment4D.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/ModelSpaceCovariance4DBase.h"
#include "oops/base/Variables.h"
#include "oops/generic/OoBump.h"
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

  std::unique_ptr<OoBump> ooBump_;
  std::vector<util::DateTime> timeslots_;
};

// =============================================================================

template<typename MODEL>
ErrorCovariance4DBUMP<MODEL>::ErrorCovariance4DBUMP(const Geometry_ & resol,
                                                    const Variables & vars,
                                                    const eckit::Configuration & conf,
                                                    const State4D_ & xb, const State4D_ & fg)
  : ModelSpaceCovariance4DBase<MODEL>(xb, fg, resol, conf), ooBump_(),
    timeslots_(xb.validTimes())
{
  Log::trace() << "ErrorCovariance4DBUMP::ErrorCovariance4DBUMP starting" << std::endl;

// Setup parameters
  Parameters_ param(resol, vars, timeslots_, conf);

// Transfer OoBump pointer
  ooBump_.reset(new OoBump(param.getOoBump()));

  Log::trace() << "ErrorCovariance4DBUMP::ErrorCovariance4DBUMP done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance4DBUMP<MODEL>::~ErrorCovariance4DBUMP() {
  Log::trace() << "ErrorCovariance4DBUMP<MODEL>::~ErrorCovariance4DBUMP starting" << std::endl;
  util::Timer timer(classname(), "~ErrorCovariance4DBUMP");
  Log::trace() << "ErrorCovariance4DBUMP<MODEL>::~ErrorCovariance4DBUMP done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance4DBUMP<MODEL>::doRandomize(Increment4D_ & dx) const {
  Log::trace() << "ErrorCovariance4DBUMP<MODEL>::doRandomize starting" << std::endl;
  util::Timer timer(classname(), "doRandomize");
  UnstructuredGrid ug(ooBump_->getColocated(), timeslots_.size());
  dx.ug_coord(ug);
  ooBump_->randomizeNicas(ug);
  dx.field_from_ug(ug);
  Log::trace() << "ErrorCovariance4DBUMP<MODEL>::doRandomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance4DBUMP<MODEL>::doMultiply(const Increment4D_ & dxi,
                                              Increment4D_ & dxo) const {
  Log::trace() << "ErrorCovariance4DBUMP<MODEL>::doMultiply starting" << std::endl;
  util::Timer timer(classname(), "doMultiply");
  UnstructuredGrid ug(ooBump_->getColocated(), timeslots_.size());
  dxi.field_to_ug(ug);
  ooBump_->multiplyNicas(ug);
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
