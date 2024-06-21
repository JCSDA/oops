/*
 * (C) Crown Copyright 2023, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSOPERATORPERT_H_
#define OOPS_BASE_OBSOPERATORPERT_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "oops/base/Locations.h"
#include "oops/base/ObsOperatorBase.h"
#include "oops/base/ObsVector.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/LinearObsOperator.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

class Variables;

// -----------------------------------------------------------------------------

/// \brief Wrapper for a LinearObsOperator under the guise of a regular ObsOperator,
/// used only in Control-Pert EDA.

template <typename OBS>
class ObsOperatorPert : public ObsOperatorBase<OBS>,
                        public util::Printable,
                        private boost::noncopyable,
                        private util::ObjectCounter<ObsOperatorPert<OBS> > {
  typedef GeoVaLs<OBS>               GeoVaLs_;
  typedef ObsDiagnostics<OBS>        ObsDiags_;
  typedef Locations<OBS>             Locations_;
  typedef LinearObsOperator<OBS>     LinObsOperator_;
  typedef ObsAuxControl<OBS>         ObsAuxControl_;
  typedef ObsAuxIncrement<OBS>       ObsAuxIncrement_;
  typedef ObsOperator<OBS>           ObsOperator_;
  typedef ObsVector<OBS>             ObsVector_;
  typedef ObsSpace<OBS>              ObsSpace_;
  typedef ObsDataVector<OBS, int>    ObsDataInt_;

 public:
  static const std::string classname() {return "oops::ObsOperatorPert";}

  /// Set up \p linoper under the guise of a regular observation operator.
  /// A regular obs operator also needs to be created with \p obsspace and \p parameters
  /// in order to access the Locations.
  ObsOperatorPert(const ObsSpace_ & obsspace, const eckit::Configuration &,
                  const LinObsOperator_ & linoper, bool trajSet = false);
  ~ObsOperatorPert();

  /// Compute forward operator \p y = LinearObsOperator (\p x).
  /// \param[in]  x        obs operator input, Increment (under the guise of a State) interpolated
  ///                      to observations locations.
  /// \param[out] y        result of computing obs operator on \p x.
  /// \param[in]  obsaux   additional input for computing H(x), e.g. bias correction coefficients
  ///                      or obs operator parameters.
  /// \params[in] qc_flags   quality control flags
  void simulateObs(const GeoVaLs_ & x_int, ObsVector_ & y, const ObsAuxControl_ & obsaux,
                   const ObsDataInt_ & qc_flags,
                   ObsVector_ &, ObsDiags_ &) const;

  /// Variables required from the model Increment (under the guise of a State) to compute
  /// the forward operator. These variables will be provided in GeoVaLs passed to simulateObs.
  const Variables & requiredVars() const;
  /// Locations used for computing GeoVaLs that will be passed to simulateObs.
  Locations_ locations() const;

  void computeReducedVars(const oops::Variables & vars, GeoVaLs_ & gvals) const {}

 private:
  /// Print, used for logging
  void print(std::ostream &) const;

  const std::string name_;
  std::unique_ptr<ObsOperator_> oper_;  // Used to access Locations and computeReducedVars
  const LinObsOperator_ & linOper_;
  bool trajectorySet;
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsOperatorPert<OBS>::ObsOperatorPert(const ObsSpace_ & obsspace, const eckit::Configuration & conf,
                                      const LinObsOperator_ & linoper, bool trajSet)
  : name_("oops::ObsOperatorPert::"+obsspace.obsname()), linOper_(linoper), trajectorySet(trajSet)
{
  Log::trace() << "ObsOperatorPert<OBS>::ObsOperatorPert starting" << std::endl;
  util::Timer timer(name_, "ObsOperatorPert");
  oper_.reset(new ObsOperator_(obsspace, conf));
  Log::trace() << "ObsOperatorPert<OBS>::ObsOperatorPert done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
ObsOperatorPert<OBS>::~ObsOperatorPert() {
  Log::trace() << "ObsOperatorPert<OBS>::~ObsOperatorPert starting" << std::endl;
  util::Timer timer(name_, "~ObsOperatorPert");
  oper_.reset();
  Log::trace() << "ObsOperatorPert<OBS>::~ObsOperatorPert done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsOperatorPert<OBS>::simulateObs(const GeoVaLs_ & gvals, ObsVector_ & yy,
                                       const ObsAuxControl_ & aux,
                                       const ObsDataInt_ & qc_flags,
                                       ObsVector_ &,
                                       ObsDiags_ &) const {
  ASSERT(trajectorySet);
  Log::trace() << "ObsOperatorPert<OBS>::simulateObs starting" << std::endl;
  util::Timer timer(name_, "simulateObs");
  ObsAuxControl_ auxcopy(aux, false);
  ObsAuxIncrement_ auxInc(auxcopy.obspace(), auxcopy.config());
  auxInc.diff(aux, auxcopy);    // In UFO implementation, auxInc does not include coefficients
                                // for static BC predictors but only variational ones
  linOper_.simulateObsTL(gvals, yy, auxInc, qc_flags);
  Log::trace() << "ObsOperatorPert<OBS>::simulateObs done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
const Variables & ObsOperatorPert<OBS>::requiredVars() const {
  Log::trace() << "ObsOperatorPert<OBS>::requiredVars starting" << std::endl;
  util::Timer timer(name_, "requiredVars");
  return linOper_.requiredVars();
}

// -----------------------------------------------------------------------------

template <typename OBS>
Locations<OBS> ObsOperatorPert<OBS>::locations() const {
  Log::trace() << "ObsOperatorPert<OBS>::locations starting" << std::endl;
  util::Timer timer(name_, "locations");
  return Locations_(oper_->locations());
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsOperatorPert<OBS>::print(std::ostream & os) const {
  Log::trace() << "ObsOperatorPert<OBS>::print starting" << std::endl;
  util::Timer timer(name_, "print");
  os << linOper_;
  Log::trace() << "ObsOperatorPert<OBS>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSOPERATORPERT_H_

