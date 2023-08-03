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

#ifndef OOPS_BASE_OBSERVERTLAD_H_
#define OOPS_BASE_OBSERVERTLAD_H_

#include <memory>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/Geometry.h"
#include "oops/base/GetValues.h"
#include "oops/base/Locations.h"
#include "oops/base/ObserverUtils.h"
#include "oops/base/ObsVector.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/LinearObsOperator.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/DateTime.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

/// Computes observation equivalent TL and AD to/from increments.
template <typename MODEL, typename OBS>
class ObserverTLAD {
  typedef Geometry<MODEL>              Geometry_;
  typedef GeoVaLs<OBS>                 GeoVaLs_;
  typedef GetValues<MODEL, OBS>        GetValues_;
  typedef LinearObsOperator<OBS>       LinearObsOperator_;
  typedef Locations<OBS>               Locations_;
  typedef ObsAuxControl<OBS>           ObsAuxCtrl_;
  typedef ObsAuxIncrement<OBS>         ObsAuxIncr_;
  typedef ObsOperator<OBS>             ObsOperator_;
  typedef ObsSpace<OBS>                ObsSpace_;
  typedef ObsVector<OBS>               ObsVector_;
  typedef ObserverParameters<OBS>      Parameters_;

 public:
  ObserverTLAD(const ObsSpace_ &, const Parameters_ &);
  ~ObserverTLAD() {}

  std::vector<std::shared_ptr<GetValues_>> initializeTraj(const Geometry_ &, const ObsAuxCtrl_ &);
  void finalizeTraj();

  void finalizeTL(const ObsAuxIncr_ &, ObsVector_ &);

  void initializeAD(const ObsVector_ &, ObsAuxIncr_ &);
  void finalizeAD() {}

/// Accessor to linear obs operator
  const LinearObsOperator_ & linObsOp() {return hoptlad_;}

 private:
  typedef std::vector<size_t> VariableSizes;

  Parameters_                   parameters_;
  const ObsSpace_ &             obspace_;    // ObsSpace used in H(x)
  Variables                     hopVars_;
  VariableSizes                 hopVarSizes_;   // Sizes of variables requested from model
  ObsOperator_                  hop_;        // Obs operator
  LinearObsOperator_            hoptlad_;    // Linear obs operator
  // Instances of GetValues. Each receives a list of model variables and a set of paths along which
  // these variables should be interpolated. The interpolated values are stored in a single GeoVaLs
  // object (shared between all instances of GetValues).
  std::vector<std::shared_ptr<GetValues_>> getvals_;
  VariableSizes        hoptladVarSizes_;     // Sizes of variables requested from model for
                                             // TL/AD (e.g. number of vertical levels)
  std::unique_ptr<Locations_> locations_;
  util::DateTime winbgn_;                    // Begining of assimilation window
  util::DateTime winend_;                    // End of assimilation window
  const ObsAuxCtrl_ *           ybias_;
  bool init_;
};

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
ObserverTLAD<MODEL, OBS>::ObserverTLAD(const ObsSpace_ & obsdb, const Parameters_ & params)
  : parameters_(params), obspace_(obsdb), hopVars_(), hopVarSizes_(),
    hop_(obspace_, parameters_.obsOperator),
    hoptlad_(obspace_,
             params.linearObsOperator.value() != boost::none ?
               params.linearObsOperator.value().value() :
               // Hack: when "linear obs operator" is not specified in the input file, reinterpret
               //       the entry for "obs operator" as a linear obs operator option. In the long
               //       term, we need a design that either,
               //       - allows constructing LinearObsOperator from either set of Parameters, or
               //       - merges the two sets of Parameters so this switch can be removed
               validateAndDeserialize<typename LinearObsOperator_::Parameters_>(
                   params.obsOperator.value().toConfiguration())),
    getvals_(), winbgn_(obsdb.windowStart()), winend_(obsdb.windowEnd()),
    ybias_(nullptr), init_(false)
{
  Log::trace() << "ObserverTLAD::ObserverTLAD" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
std::vector<std::shared_ptr<GetValues<MODEL, OBS>>>
ObserverTLAD<MODEL, OBS>::initializeTraj(const Geometry_ & geom, const ObsAuxCtrl_ & ybias) {
  Log::trace() << "ObserverTLAD::initializeTraj start" << std::endl;
  ybias_ = &ybias;

  // Get the list of variables to be obtained from the model state
  oops::Variables geovars = hop_.requiredVars();
  geovars += ybias_->requiredVars();

  // Get the observation locations and their discretizations
  locations_ = std::make_unique<Locations_>(hop_.locations());

  // Variables grouped by the set of paths along which they'll be interpolated
  const std::vector<oops::Variables> groupedHopVars = groupVariablesByLocationSamplingMethod(
        geovars, *locations_);

  // Get the sizes of these variables (i.e. the numbers of levels in the corresponding GeoVaLs)
  const std::vector<VariableSizes> groupedHopVarSizes = variableSizes(groupedHopVars, geom);

  std::tie(hopVars_, hopVarSizes_) = mergeVariablesAndSizes(groupedHopVars, groupedHopVarSizes);

  // Now deal with the variables obtained from the model increment.
  hoptladVarSizes_ = geom.variableSizes(hoptlad_.requiredVars());
  const std::vector<Variables> groupedHoptladVars = groupVariablesByLocationSamplingMethod(
        hoptlad_.requiredVars(), *locations_);

  // Set up GetValues
  getvals_ = makeGetValuesVector(parameters_.getValues.value(), geom, winbgn_, winend_,
                                 *locations_, groupedHopVars, groupedHoptladVars);

  init_ = true;
  Log::trace() << "ObserverTLAD::initializeTraj done" << std::endl;
  return getvals_;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserverTLAD<MODEL, OBS>::finalizeTraj() {
  Log::trace() << "ObserverTLAD::finalizeTraj start" << std::endl;
  ASSERT(init_);

  // Fill geovals
  GeoVaLs_ geovals(*locations_, hopVars_, hopVarSizes_);
  for (size_t m = 0; m < getvals_.size(); ++m) {
    getvals_[m]->fillGeoVaLs(geovals);
  }

  // Compute the reduced representation of the GeoVaLs for which it's been requested
  oops::Variables reducedVars = ybias_->requiredVars();
  hop_.computeReducedVars(reducedVars, geovals);

  /// Set linearization trajectory for H(x)
  hoptlad_.setTrajectory(geovals, *ybias_);

  init_ = false;
  Log::trace() << "ObserverTLAD::finalizeTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserverTLAD<MODEL, OBS>::finalizeTL(const ObsAuxIncr_ & ybiastl, ObsVector_ & ydeptl) {
  Log::trace() << "ObserverTLAD::finalizeTL start" << std::endl;

  // TODO(wsmigaj): should we allow linear operators to require also *reduced* GeoVaLs?
  GeoVaLs_ geovals(*locations_, hoptlad_.requiredVars(), hoptladVarSizes_);

  // Fill GeoVaLs
  for (size_t m = 0; m < getvals_.size(); ++m) {
    getvals_[m]->fillGeoVaLsTL(geovals);
  }

  // Compute linear H(x)
  hoptlad_.simulateObsTL(geovals, ydeptl, ybiastl);

  Log::trace() << "ObserverTLAD::finalizeTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserverTLAD<MODEL, OBS>::initializeAD(const ObsVector_ & ydepad, ObsAuxIncr_ & ybiasad) {
  Log::trace() << "ObserverTLAD::initializeAD start" << std::endl;

  GeoVaLs_ geovals(*locations_, hoptlad_.requiredVars(), hoptladVarSizes_);

  // Compute adjoint of H(x)
  hoptlad_.simulateObsAD(geovals, ydepad, ybiasad);

  // GeoVaLs forcing to GetValues
  for (size_t m = 0; m < getvals_.size(); ++m) {
    getvals_[m]->fillGeoVaLsAD(geovals);
  }

  Log::trace() << "ObserverTLAD::initializeAD done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVERTLAD_H_
