/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_ERRORCOVARIANCE_H_
#define OOPS_INTERFACE_ERRORCOVARIANCE_H_

#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "eckit/system/ResourceUsage.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/Increment4D.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

// Should factory be here and generic covariances wrtten at the MODEL::Increment level? YT

/// Wrapper for model space error covariances.
/*!
 *  This class provides the operations associated with the model space error
 *  covariance matrices (B or Q). It wraps the actual error covariance matrix
 *  which can be a model specific one or a generic one.
 */

template <typename MODEL>
class ErrorCovariance : public oops::ModelSpaceCovarianceBase<MODEL>,
                        public util::Printable,
                        private util::ObjectCounter<ErrorCovariance<MODEL> >,
                        private boost::noncopyable {
  typedef typename MODEL::Covariance Covariance_;
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment4D<MODEL>         Increment4D_;
  typedef Increment<MODEL>           Increment_;
  typedef State4D<MODEL>             State_;

 public:
  /// Defined as Covariance_::Parameters_ if Covariance_ defines a Parameters_ type; otherwise as
  /// GenericModelSpaceCovarianceParameters<MODEL>.
  typedef TParameters_IfAvailableElseFallbackType_t<
    Covariance_, GenericModelSpaceCovarianceParameters<MODEL>> Parameters_;

  static const std::string classname() {return "oops::ErrorCovariance";}

  ErrorCovariance(const Geometry_ &, const Variables &, const Parameters_ &,
                  const State_ &, const State_ &);
  ErrorCovariance(const Geometry_ &, const Variables &, const eckit::Configuration &,
                  const State_ &, const State_ &);
  virtual ~ErrorCovariance();

 private:
  void doRandomize(Increment4D_ &) const override;
  void doMultiply(const Increment4D_ &, Increment4D_ &) const override;
  void doInverseMultiply(const Increment4D_ &, Increment4D_ &) const override;

  void print(std::ostream &) const override;

  std::unique_ptr<Covariance_> covariance_;
};

// =============================================================================

template<typename MODEL>
ErrorCovariance<MODEL>::ErrorCovariance(const Geometry_ & resol, const Variables & vars,
                                        const Parameters_ & parameters,
                                        const State_ & xb, const State_ & fg)
  : ModelSpaceCovarianceBase<MODEL>(resol, parameters, xb, fg),
    covariance_()
{
  Log::trace() << "ErrorCovariance<MODEL>::ErrorCovariance starting" << std::endl;
  util::Timer timer(classname(), "ErrorCovariance");
  size_t init = eckit::system::ResourceUsage().maxResidentSetSize();
  covariance_.reset(new Covariance_(resol.geometry(), vars,
                                    parametersOrConfiguration<HasParameters_<Covariance_>::value>(
                                      parameters),
                                    xb[0].state(), fg[0].state()));
  size_t current = eckit::system::ResourceUsage().maxResidentSetSize();
  this->setObjectSize(current - init);
  Log::trace() << "ErrorCovariance<MODEL>::ErrorCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance<MODEL>::ErrorCovariance(const Geometry_ & resol, const Variables & vars,
                                        const eckit::Configuration & conf,
                                        const State_ & xb, const State_ & fg)
  : ErrorCovariance<MODEL>(resol, vars, conf, xb, fg)
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance<MODEL>::~ErrorCovariance() {
  Log::trace() << "ErrorCovariance<MODEL>::~ErrorCovariance starting" << std::endl;
  util::Timer timer(classname(), "~ErrorCovariance");
  covariance_.reset();
  Log::trace() << "ErrorCovariance<MODEL>::~ErrorCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::doRandomize(Increment4D_ & dx) const {
  Log::trace() << "ErrorCovariance<MODEL>::doRandomize starting" << std::endl;
  util::Timer timer(classname(), "doRandomize");
  covariance_->randomize(dx[0].increment());
  Log::trace() << "ErrorCovariance<MODEL>::doRandomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::doMultiply(const Increment4D_ & dx_in, Increment4D_ & dx_out) const {
  Log::trace() << "ErrorCovariance<MODEL>::doMultiply starting" << std::endl;
  util::Timer timer(classname(), "doMultiply");
  const std::vector<util::DateTime> times(dx_in.validTimes());
  Increment4D_ tmpInc = dx_in;
  for (size_t sw = 1; sw < tmpInc.local_time_size(); sw++) {
    tmpInc[0].axpy(1.0, tmpInc[sw], false);
  }
  if (tmpInc.commTime().rank() > 0) {
    oops::mpi::send(tmpInc.commTime(), tmpInc[0], 0, 0);
  } else {
    Increment_ dxtmp(tmpInc[0], false);
    for (size_t sw = 1; sw < dx_out.commTime().size(); ++sw) {
      // collect all the sums (dxtmp) from each rank
      oops::mpi::receive(tmpInc.commTime(), dxtmp, sw, 0);
      tmpInc[0].axpy(1.0, dxtmp, false);
    }
    covariance_->multiply(tmpInc[0].increment(), dx_out[0].increment());
  }
  // Broadcast sum to other ranks
  oops::mpi::broadcast(tmpInc.commTime(), dx_out[0], 0);
  // Copy increment to all sub-windows
  dx_out[0].updateTime(times[0] - dx_out[0].validTime());
  for (size_t jt = 1; jt < dx_out.local_time_size(); ++jt) {
    dx_out[jt] = dx_out[0];
    dx_out[jt].updateTime(times[jt] - dx_out[jt].validTime());
  }
  Log::trace() << "ErrorCovariance<MODEL>::doMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::doInverseMultiply(const Increment4D_ & dx1, Increment4D_ & dx2) const {
  Log::trace() << "ErrorCovariance<MODEL>::doInverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "doInverseMultiply");
  covariance_->inverseMultiply(dx1[0].increment(), dx2[0].increment());
  Log::trace() << "ErrorCovariance<MODEL>::doInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ErrorCovariance<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *covariance_;
  Log::trace() << "ErrorCovariance<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_ERRORCOVARIANCE_H_
