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

#include "oops/base/ObsVector.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

class Variables;

// -----------------------------------------------------------------------------
/// \brief MODEL-agnostic part of nonlinear observation (forward) operator.
/// The full nonlinear observation operator from State x to ObsVector is:
/// ObsOperator ( GetValues (State) )
/// ObsOperator uses GeoVaLs (result of GetValues(State) - model State at
/// observations locations) as input data to compute forward operator.
///
/// Note: each implementation should typedef `Parameters_` to the name of a subclass of
/// oops::Parameters holding its configuration settings and provide a constructor with the
/// following signature:
///
///     ObsOperator(const OBS::ObsSpace &, const Parameters_ &);
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
  /// A subclass of oops::Parameters holding the configuration settings of the operator.
  typedef typename ObsOperator_::Parameters_ Parameters_;

  static const std::string classname() {return "oops::ObsOperator";}

  /// Set up observation operator for the \p obsspace observations, with
  /// parameters defined in \p parameters
  ObsOperator(const ObsSpace_ & obsspace, const Parameters_ & parameters);
  ~ObsOperator();

  /// Compute forward operator \p y = ObsOperator (\p x).
  /// \param[in]  x        obs operator input, State interpolated to observations locations.
  /// \param[out] y        result of computing obs operator on \p x.
  /// \param[in]  obsaux   additional input for computing H(x), used in the minimization
  ///                      in Variational DA, e.g. bias correction coefficients or obs operator
  ///                      parameters.
  /// \param[out] obsbias  bias correction of the departure between \p y and the observed values;
  ///                      when \p obsbias is non-zero, it is added to \p y within the obs
  ///                      operator
  /// \param[out] obsdiags   additional diagnostics output from computing obs operator that is not
  ///                        used in the assimilation, and can be used by ObsFilters.
  void simulateObs(const GeoVaLs_ & x_int, ObsVector_ & y, const ObsAuxControl_ & obsaux,
                   ObsVector_ & obsbias, ObsDiags_ & obsdiags) const;

  /// Variables required from the model State to compute obs operator. These variables
  /// will be provided in GeoVaLs passed to simulateObs.
  const Variables & requiredVars() const;
  /// Locations used for computing GeoVaLs that will be passed to simulateObs.
  Locations_ locations() const;

 private:
  /// Print, used for logging
  void print(std::ostream &) const;

  const std::string name_;
  /// Pointer to the implementation of ObsOperator
  std::unique_ptr<ObsOperator_> oper_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsOperator<OBS>::ObsOperator(const ObsSpace_ & os, const Parameters_ & parameters)
  : name_("oops::ObsOperator::"+os.obsname()), oper_()
{
  Log::trace() << "ObsOperator<OBS>::ObsOperator starting" << std::endl;
  util::Timer timer(name_, "ObsOperator");
  oper_.reset(new ObsOperator_(os.obsspace(), parameters));
  Log::trace() << "ObsOperator<OBS>::ObsOperator done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
ObsOperator<OBS>::~ObsOperator() {
  Log::trace() << "ObsOperator<OBS>::~ObsOperator starting" << std::endl;
  util::Timer timer(name_, "~ObsOperator");
  oper_.reset();
  Log::trace() << "ObsOperator<OBS>::~ObsOperator done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsOperator<OBS>::simulateObs(const GeoVaLs_ & gvals, ObsVector_ & yy,
                                     const ObsAuxControl_ & aux, ObsVector_ & ybias,
                                     ObsDiags_ & ydiag) const {
  Log::trace() << "ObsOperator<OBS>::simulateObs starting" << std::endl;
  util::Timer timer(name_, "simulateObs");
  oper_->simulateObs(gvals.geovals(), yy.obsvector(), aux.obsauxcontrol(), ybias.obsvector(),
                     ydiag.obsdiagnostics());
  Log::trace() << "ObsOperator<OBS>::simulateObs done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
const Variables & ObsOperator<OBS>::requiredVars() const {
  Log::trace() << "ObsOperator<OBS>::requiredVars starting" << std::endl;
  util::Timer timer(name_, "requiredVars");
  return oper_->requiredVars();
}

// -----------------------------------------------------------------------------

template <typename OBS>
Locations<OBS> ObsOperator<OBS>::locations() const {
  Log::trace() << "ObsOperator<OBS>::locations starting" << std::endl;
  util::Timer timer(name_, "locations");
  return Locations_(oper_->locations());
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsOperator<OBS>::print(std::ostream & os) const {
  Log::trace() << "ObsOperator<OBS>::print starting" << std::endl;
  util::Timer timer(name_, "print");
  os << *oper_;
  Log::trace() << "ObsOperator<OBS>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSOPERATOR_H_
