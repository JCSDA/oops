/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_LINEAROBSOPERATOR_H_
#define OOPS_INTERFACE_LINEAROBSOPERATOR_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "oops/base/ObsVector.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

class Variables;

// -----------------------------------------------------------------------------
/// \brief MODEL-agnostic part of tangent-linear and adjoint of the nonlinear
/// observation (forward) operator ObsOperator.
///
/// Note: each implementation should typedef `Parameters_` to the name of a subclass of
/// oops::Parameters holding its configuration settings and provide a constructor with the
/// following signature:
///
///     LinearObsOperator(const OBS::ObsSpace &, const Parameters_ &);
template <typename OBS>
class LinearObsOperator : public util::Printable,
                          private boost::noncopyable,
                          private util::ObjectCounter<LinearObsOperator<OBS> > {
  typedef typename OBS::LinearObsOperator  LinearObsOper_;
  typedef GeoVaLs<OBS>             GeoVaLs_;
  typedef ObsAuxControl<OBS>       ObsAuxControl_;
  typedef ObsAuxIncrement<OBS>     ObsAuxIncrement_;
  typedef ObsSpace<OBS>            ObsSpace_;
  typedef ObsVector<OBS>           ObsVector_;
  typedef ObsDataVector<OBS, int>  ObsDataInt_;

 public:
  static const std::string classname() {return "oops::LinearObsOperator";}

  /// Set up TL and AD of observation operator for the \p obsspace observations, with
  /// parameters defined in \p parameters.
  LinearObsOperator(const ObsSpace_ & obsspace, const eckit::Configuration &);
  ~LinearObsOperator();

  /// Sets up the trajectory for future calls of simulateObsTL or simulateObsAD.
  /// The implementations could e.g. save the trajectory \p x0, or compute and save the Jacobian
  /// of observation operator around \p x0.
  /// Always called before simulateObsTL or simulateObsAD.
  /// \param[in]  x0       trajectory for linearization of obs operator, State interpolated
  ///                      to observations locations (defined by ObsOperator::locations())
  /// \param[in]  obsaux   additional obs operator input, used in the minimization
  ///                      in Variational DA, e.g. bias correction coefficients or obs operator
  ///                      parameters.
  void setTrajectory(const GeoVaLs_ & x0, const ObsAuxControl_ & obsaux,
                     const ObsDataInt_ & qc_flags);

  /// Apply tangent-linear of the observation operator linearized around the trajectory that was
  /// passed to setTrajectory method (which is always called before simulateObsTL).
  /// \param[in]  dx       input to the TL obs operator, Increment interpolated to observations
  ///                      locations.
  /// \param[out] dy       output of the TL obs operator.
  /// \param[in]  dobsaux: additional input to the TL obs operator, e.g. perturbation to bias
  ///                      coefficients or obs operator parameters
  /// \param[in] qc_flags quality control flags
  void simulateObsTL(const GeoVaLs_ & dx, ObsVector_ & dy, const ObsAuxIncrement_ & dobsaux,
                     const ObsDataInt_ & qc_flags) const;
  /// Apply adjoint of the observation operator linearized around the trajectory that was
  /// passed to setTrajectory method (which is always called before simulateObsAD).
  /// \param[out] dx       output of the AD obs operator, Increment interpolated to observations
  ///                      locations.
  /// \param[in]  dy       input of the AD obs operator, perturbation to the ObsVector.
  /// \param[out] dobsaux  additional output of the AD obs operator, e.g. perturbation to bias
  ///                      coefficients or obs operator parameters.
  /// \param[in] qc_flags quality control flags
  void simulateObsAD(GeoVaLs_ & dx, const ObsVector_ & dy, ObsAuxIncrement_ & dobsaux,
                     const ObsDataInt_ & qc_flags) const;
  /// Variables required from the model Increment to compute TL or AD of the obs operator.
  /// These variables will be provided in GeoVaLs passed to simulateObsTL and simulateObsAD.
  /// Note: these Variables may be different from variables returned by ObsOperator::requiredVars(),
  /// which will be provided in GeoVaLs passed to setTrajectory.
  const Variables & requiredVars() const;

 private:
  /// Print, used for logging
  void print(std::ostream &) const;

  const std::string name_;
  /// Pointer to the implementation of LinearObsOperator
  std::unique_ptr<LinearObsOper_> oper_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
LinearObsOperator<OBS>::LinearObsOperator(const ObsSpace_ & os, const eckit::Configuration & config)
  : name_("oops::LinearObsOper::"+os.obsname()), oper_()
{
  Log::trace() << "LinearObsOperator<OBS>::LinearObsOperator starting" << std::endl;
  util::Timer timer(name_, "LinearObsOperator");
  oper_.reset(new LinearObsOper_(os.obsspace(), config));
  Log::trace() << "LinearObsOperator<OBS>::LinearObsOperator done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
LinearObsOperator<OBS>::~LinearObsOperator() {
  Log::trace() << "LinearObsOperator<OBS>::~LinearObsOperator starting" << std::endl;
  util::Timer timer(name_, "~LinearObsOperator");
  oper_.reset();
  Log::trace() << "LinearObsOperator<OBS>::~LinearObsOperator done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void LinearObsOperator<OBS>::setTrajectory(const GeoVaLs_ & gvals, const ObsAuxControl_ & aux,
                                           const ObsDataInt_ & qc_flags) {
  Log::trace() << "LinearObsOperator<OBS>::setTrajectory starting" << std::endl;
  util::Timer timer(name_, "setTrajectory");
  oper_->setTrajectory(gvals.geovals(), aux.obsauxcontrol(),
                       qc_flags.obsdatavector());
  Log::trace() << "LinearObsOperator<OBS>::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void LinearObsOperator<OBS>::simulateObsTL(const GeoVaLs_ & gvals, ObsVector_ & yy,
                                             const ObsAuxIncrement_ & aux,
                                             const ObsDataInt_ & qc_flags) const {
  Log::trace() << "LinearObsOperator<OBS>::simulateObsTL starting" << std::endl;
  util::Timer timer(name_, "simulateObsTL");
  oper_->simulateObsTL(gvals.geovals(), yy.obsvector(), aux.obsauxincrement(),
                       qc_flags.obsdatavector());
  Log::trace() << "LinearObsOperator<OBS>::simulateObsTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void LinearObsOperator<OBS>::simulateObsAD(GeoVaLs_ & gvals, const ObsVector_ & yy,
                                             ObsAuxIncrement_ & aux,
                                             const ObsDataInt_ & qc_flags) const {
  Log::trace() << "LinearObsOperator<OBS>::simulateObsAD starting" << std::endl;
  util::Timer timer(name_, "simulateObsAD");
  oper_->simulateObsAD(gvals.geovals(), yy.obsvector(), aux.obsauxincrement(),
                       qc_flags.obsdatavector());
  Log::trace() << "LinearObsOperator<OBS>::simulateObsAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
const Variables & LinearObsOperator<OBS>::requiredVars() const {
  Log::trace() << "LinearObsOperator<OBS>::requiredVars starting" << std::endl;
  util::Timer timer(name_, "requiredVars");
  return oper_->requiredVars();
}

// -----------------------------------------------------------------------------

template<typename OBS>
void LinearObsOperator<OBS>::print(std::ostream & os) const {
  Log::trace() << "LinearObsOperator<OBS>::print starting" << std::endl;
  util::Timer timer(name_, "print");
  os << *oper_;
  Log::trace() << "LinearObsOperator<OBS>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_LINEAROBSOPERATOR_H_
