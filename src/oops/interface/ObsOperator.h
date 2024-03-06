/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Crown Copyright 2023, the Met Office.
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

#include "oops/base/Locations.h"
#include "oops/base/ObsOperatorBase.h"
#include "oops/base/ObsVector.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsDataVector.h"
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
class ObsOperator : public ObsOperatorBase<OBS>,
                    public util::Printable,
                    private boost::noncopyable,
                    private util::ObjectCounter<ObsOperator<OBS> > {
  typedef typename OBS::ObsOperator  ObsOperator_;
  typedef GeoVaLs<OBS>               GeoVaLs_;
  typedef Locations<OBS>             Locations_;
  typedef ObsDiagnostics<OBS>        ObsDiags_;
  typedef ObsAuxControl<OBS>         ObsAuxControl_;
  typedef ObsVector<OBS>             ObsVector_;
  typedef ObsSpace<OBS>              ObsSpace_;
  typedef ObsDataVector<OBS, int>    ObsDataInt_;

 public:
  static const std::string classname() {return "oops::ObsOperator";}

  /// Set up observation operator for the \p obsspace observations, with
  /// parameters defined in \p parameters
  ObsOperator(const ObsSpace_ & obsspace, const eckit::Configuration &);
  ~ObsOperator();

  /// Compute forward operator \p y = ObsOperator (\p x).
  /// \param[in]  x        obs operator input, State interpolated along paths sampling the
  ///                      observation locations.
  /// \param[out] y        result of computing obs operator on \p x.
  /// \param[in]  obsaux   additional input for computing H(x), used in the minimization
  ///                      in Variational DA, e.g. bias correction coefficients or obs operator
  ///                      parameters.
  /// \param[out] obsbias  bias correction of the departure between \p y and the observed values;
  ///                      when \p obsbias is non-zero, it is added to \p y within the obs
  ///                      operator
  /// \param[out] obsdiags   additional diagnostics output from computing obs operator that is not
  ///                        used in the assimilation, and can be used by ObsFilters.
  /// \param[in] qc_flags  quality control flags
  void simulateObs(const GeoVaLs_ & x_int, ObsVector_ & y, const ObsAuxControl_ & obsaux,
                   const ObsDataInt_ & qc_flags,
                   ObsVector_ & obsbias, ObsDiags_ & obsdiags) const;

  /// Variables required from the model State to compute the obs operator.
  const Variables & requiredVars() const;

  /// Return an object holding one or more collections of paths sampling the observation locations
  /// and indicating along which of these sets of paths individual model variables should be
  /// interpolated.
  ///
  /// Operators simulating pointwise observations will normally ask for all variables to be
  /// interpolated along the same set of vertical paths, each at the nominal latitude
  /// and longitude of a single observation. Operators simulating spatially extended observations
  /// can ask for some or all variables to be interpolated along multiple paths sampling each
  /// observation location, for instance to then compute a weighted average.
  Locations_ locations() const;

  /// \brief Convert values of model variables stored in the sampled format to the reduced format.
  ///
  /// This typically consists in computing, for each location, a weighted average of all the
  /// profiles obtained by model field interpolation along the paths sampling that location.
  ///
  /// \param[in] vars
  ///   List of variables whose reduced representation should be computed.
  /// \param[inout] gvals
  ///   A container for the sampled and reduced representations of the values of model variables.
  ///   Values stored in the sampled format will already have been filled in by sampling model
  ///   fields along the interpolation paths specified by locations(). This function needs to fill
  ///   in the reduced representation of at least the variables `vars` (unless it is already
  ///   available -- some implementations of the GeoVaLs interface automatically detect variables
  ///   whose sampled and reduced formats are identical, store only their sampled representation
  ///   and treat the reduced representation as an alias for the sampled one). It may optionally
  ///   compute the reduced representation of other variables as well.
  void computeReducedVars(const oops::Variables & vars, GeoVaLs_ & gvals) const;

 private:
  /// Print, used for logging
  void print(std::ostream &) const;

  const std::string name_;
  /// Pointer to the implementation of ObsOperator
  std::unique_ptr<ObsOperator_> oper_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsOperator<OBS>::ObsOperator(const ObsSpace_ & os, const eckit::Configuration & config)
  : name_("oops::ObsOperator::"+os.obsname()), oper_()
{
  Log::trace() << "ObsOperator<OBS>::ObsOperator starting" << std::endl;
  util::Timer timer(name_, "ObsOperator");
  oper_.reset(new ObsOperator_(os.obsspace(), config));
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
                                     const ObsAuxControl_ & aux,
                                     const ObsDataInt_ & qc_flags,
                                     ObsVector_ & ybias,
                                     ObsDiags_ & ydiag) const {
  Log::trace() << "ObsOperator<OBS>::simulateObs starting" << std::endl;
  util::Timer timer(name_, "simulateObs");
  oper_->simulateObs(gvals.geovals(), yy.obsvector(), aux.obsauxcontrol(),
                     qc_flags.obsdatavector(),
                     ybias.obsvector(), ydiag.obsdiagnostics());
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
  return oper_->locations();
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsOperator<OBS>::computeReducedVars(const oops::Variables & vars, GeoVaLs_ & gvals) const {
  Log::trace() << "ObsOperator<OBS>::computeReducedVars starting" << std::endl;
  util::Timer timer(name_, "computeReducedVars");
  oper_->computeReducedVars(vars, gvals.geovals());
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
