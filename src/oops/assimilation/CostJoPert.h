/*
 * (C) Crown Copyright 2023, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_COSTJOPERT_H_
#define OOPS_ASSIMILATION_COSTJOPERT_H_

#include <memory>
#include <utility>
#include <vector>

#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostJo.h"
#include "oops/base/Geometry.h"
#include "oops/base/Observations.h"
#include "oops/base/Observers.h"
#include "oops/base/ObserversTLAD.h"
#include "oops/base/ObsOperatorBase.h"
#include "oops/base/ObsOperatorPert.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Jo Cost Function for the Pert members of Control-Pert EDA
/*!
 * The CostJoPert class encapsulates the Jo term of the cost function
 * in the case of the Pert members of Control-Pert EDA. The Observers
 * to be called during the model integration is managed inside the
 * CostJoPert class.
 */

template<typename MODEL, typename OBS> class CostJoPert : public CostJo<MODEL, OBS> {
  typedef ControlVariable<MODEL, OBS>   CtrlVar_;
  typedef Geometry<MODEL>               Geometry_;
  typedef Observations<OBS>             Observations_;
  typedef Observers<MODEL, OBS>         Observers_;
  typedef ObserversTLAD<MODEL, OBS>     ObserversTLAD_;
  typedef ObsOperatorBase<OBS>          ObsOperatorBase_;
  typedef ObsOperatorPert<OBS>          ObsOperatorPert_;
  typedef PostProcessorTLAD<MODEL>      PostProcTLAD_;

 public:
  /// Construct \f$ J_o\f$ from \f$ R\f$ and \f$ y_{obs}\f$.
  CostJoPert(const eckit::Configuration &, const eckit::mpi::Comm &,
             const util::DateTime &, const util::DateTime &,
             const std::shared_ptr<ObserversTLAD_> &,
             const eckit::mpi::Comm & ctime = oops::mpi::myself());

  /// Destructor
  virtual ~CostJoPert();

  /// setPostProcTraj and computeCostTraj are empty methods
  void setPostProcTraj(const CtrlVar_ &, const eckit::Configuration &,
                       const Geometry_ &, PostProcTLAD_ &) override {};
  void computeCostTraj() override {};

 private:
  /// Create zero-valued observations which will later be perturbed
  void newObs() override;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostJoPert<MODEL, OBS>::CostJoPert(const eckit::Configuration & joConf,
                                   const eckit::mpi::Comm & comm, const util::DateTime & winbgn,
                                   const util::DateTime & winend,
                                   const std::shared_ptr<ObserversTLAD_> & linObsTLAD,
                                   const eckit::mpi::Comm & ctime)
  : CostJo<MODEL, OBS>::CostJo(comm, joConf, winbgn, winend, ctime)
{
  ObserversParameters<MODEL, OBS> joparams{};
  joparams.deserialize(joConf);
//  Reset the observation operators to be the linear observation operators contained in *linObsTLAD,
//  wrapped inside ObsOperatorPert
  std::vector<std::unique_ptr<ObsOperatorBase_>> obsOpBases_;
  for (std::size_t ii = 0; ii < this->obspaces().size(); ++ii) {
    obsOpBases_.push_back(std::make_unique<ObsOperatorPert_>(this->obspaces()[ii],
                                    observerParameters(joparams.observers.value())[ii].obsOperator,
                                    (*linObsTLAD)[ii].linObsOp(), true));
  }
  this->observers().reset(new Observers_(this->obspaces(), joConf, std::move(obsOpBases_)));
//  Set the linear observers to be *linObsTLAD
  this->observersTLAD() = linObsTLAD;
  Log::trace() << "CostJoPert::CostJoPert" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
CostJoPert<MODEL, OBS>::~CostJoPert() {
  Log::trace() << "CostJoPert::~CostJoPert" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJoPert<MODEL, OBS>::newObs() {
  this->yobs().reset(new Observations_(this->obspaces(), ""));
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJOPERT_H_
