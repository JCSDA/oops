/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_OBSOPERATOR_H_
#define OOPS_INTERFACE_OBSOPERATOR_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \brief MODEL-agnostic part of nonlinear observation (forward) operator.
/// The full nonlinear observation operator from State x to ObsVector is:
/// ObsOperator ( GetValues (State) )
/// ObsOperator uses GeoVaLs (result of GetValues(State) - model State at
/// observations locations) as input data to compute forward operator.
template <typename OBS>
class ObsOperator : public util::Printable,
                    private boost::noncopyable,
                    private util::ObjectCounter<ObsOperator<OBS> > {
  typedef typename OBS::ObsOperator  ObsOperator_;
  typedef GeoVaLs<OBS>               GeoVaLs_;
  typedef ObsDiagnostics<OBS>        ObsDiags_;
  typedef Locations<OBS>             Locations_;
  typedef ObsAuxControl<OBS>         ObsAuxControl_;
  typedef ObsVector<OBS>             ObsVector_;
  typedef ObsSpace<OBS>              ObsSpace_;

 public:
  static const std::string classname() {return "oops::ObsOperator";}

  /// Set up observation operator for the \p obsspace observations, with
  /// parameters defined in \p config
  ObsOperator(const ObsSpace_ & obsspace, const eckit::Configuration & config);
  ~ObsOperator();

  /// Compute forward operator \p y = ObsOperator (\p x).
  /// \param[in]  x        obs operator input, State interpolated to observations locations.
  /// \param[out] y        result of computing obs operator on \p x.
  /// \param[in]  obsaux   additional input for computing H(x), used in the minimization
  ///                      in Variational DA, e.g. bias correction coefficients or obs operator
  ///                      parameters.
  /// \param[out] obsdiags   additional diagnostics output from computing obs operator that is not
  ///                        used in the assimilation, and can be used by ObsFilters.
  void simulateObs(const GeoVaLs_ & x_int, ObsVector_ & y, const ObsAuxControl_ & obsaux,
                   ObsDiags_ & obsdiags) const;

  /// Variables required from the model State to compute obs operator. These variables
  /// will be provided in GeoVaLs passed to simulateObs.
  const Variables & requiredVars() const;
  /// Locations used for computing GeoVaLs that will be passed to simulateObs.
  Locations_ locations() const;

 private:
  /// Print, used for logging
  void print(std::ostream &) const;

  /// Pointer to the implementation of ObsOperator
  std::unique_ptr<ObsOperator_> oper_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsOperator<OBS>::ObsOperator(const ObsSpace_ & os,
                                const eckit::Configuration & config) : oper_() {
  Log::trace() << "ObsOperator<OBS>::ObsOperator starting" << std::endl;
  util::Timer timer(classname(), "ObsOperator");
  oper_.reset(new ObsOperator_(os.obsspace(), config));
  Log::trace() << "ObsOperator<OBS>::ObsOperator done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
ObsOperator<OBS>::~ObsOperator() {
  Log::trace() << "ObsOperator<OBS>::~ObsOperator starting" << std::endl;
  util::Timer timer(classname(), "~ObsOperator");
  oper_.reset();
  Log::trace() << "ObsOperator<OBS>::~ObsOperator done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsOperator<OBS>::simulateObs(const GeoVaLs_ & gvals, ObsVector_ & yy,
                                     const ObsAuxControl_ & aux, ObsDiags_ & ydiag) const {
  Log::trace() << "ObsOperator<OBS>::simulateObs starting" << std::endl;
  util::Timer timer(classname(), "simulateObs");
  oper_->simulateObs(gvals.geovals(), yy.obsvector(), aux.obsauxcontrol(), ydiag.obsdiagnostics());
  Log::trace() << "ObsOperator<OBS>::simulateObs done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
const Variables & ObsOperator<OBS>::requiredVars() const {
  Log::trace() << "ObsOperator<OBS>::requiredVars starting" << std::endl;
  util::Timer timer(classname(), "requiredVars");
  return oper_->requiredVars();
}

// -----------------------------------------------------------------------------

template <typename OBS>
Locations<OBS> ObsOperator<OBS>::locations() const {
  Log::trace() << "ObsOperator<OBS>::locations starting" << std::endl;
  util::Timer timer(classname(), "locations");
  return Locations_(oper_->locations());
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsOperator<OBS>::print(std::ostream & os) const {
  Log::trace() << "ObsOperator<OBS>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *oper_;
  Log::trace() << "ObsOperator<OBS>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSOPERATOR_H_
