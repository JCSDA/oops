/*
 * (C) Crown Copyright 2023, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSOPERATORBASE_H_
#define OOPS_BASE_OBSOPERATORBASE_H_

#include "oops/base/Locations.h"
#include "oops/base/ObsVector.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/LinearObsOperator.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsSpace.h"

namespace oops {

class Variables;

// -----------------------------------------------------------------------------

/// Abstract base class for observation operators
/*!
 * Derived classes include oops::interface::ObsOperator and oops::ObsOperatorPert.
 * The latter is used in Control-Pert EDA.
 */

template <typename OBS>
class ObsOperatorBase {
  typedef GeoVaLs<OBS>               GeoVaLs_;
  typedef LinearObsOperator<OBS>     LinearObsOperator_;
  typedef ObsDiagnostics<OBS>        ObsDiags_;
  typedef Locations<OBS>             Locations_;
  typedef ObsAuxControl<OBS>         ObsAuxControl_;
  typedef ObsVector<OBS>             ObsVector_;
  typedef ObsSpace<OBS>              ObsSpace_;
  typedef ObsDataVector<OBS, int>    ObsDataInt_;

 public:
  virtual ~ObsOperatorBase() {}

  virtual void simulateObs(const GeoVaLs_ & x_int, ObsVector_ & y, const ObsAuxControl_ & obsaux,
                           const ObsDataInt_ & qc_flags,
                           ObsVector_ & obsbias, ObsDiags_ & obsdiags) const = 0;

  /// Variables required from the model State to compute obs operator. These variables
  /// will be provided in GeoVaLs passed to simulateObs.
  virtual const Variables & requiredVars() const = 0;

  /// Locations used for computing GeoVaLs that will be passed to simulateObs.
  virtual Locations_ locations() const = 0;

  virtual void computeReducedVars(const oops::Variables & vars, GeoVaLs_ & gvals) const = 0;
};

// -----------------------------------------------------------------------------

};  // namespace oops

#endif  // OOPS_BASE_OBSOPERATORBASE_H_
