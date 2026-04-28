/*
 * (C) Copyright 2026-2026 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_SEQUENTIALENSEMBLESOLVER_H_
#define OOPS_ASSIMILATION_SEQUENTIALENSEMBLESOLVER_H_

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/Geometry.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observers.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/StateSet.h"

#include "oops/assimilation/LocalEnsembleSolver.h"

#include "oops/util/abor1_cpp.h"

namespace oops {


/// This implements the parallel algorithm for sequential ensemble DA as
/// described in:
///
/// Anderson, J. L., & Collins, N. (2007). Scalable implementations of ensemble
/// filter algorithms for data assimilation. Journal of Atmospheric and Oceanic
/// Technology, 24(8), 1452-1463.
///
/// TODO(Travis) We are using LocalEnsembleSolver as base class for SequentialEnsembleSolver
/// mainly to  use the hofx code. Instead, LocalEnsembleSolver should have common
/// components (which is most of it), pulled out into an EnsembleSolver class
template <typename MODEL, typename OBS>
class SequentialEnsembleSolver : public LocalEnsembleSolver<MODEL, OBS> {
  typedef Departures<OBS>             Departures_;
  typedef Geometry<MODEL>             Geometry_;
  typedef GeometryIterator<MODEL>     GeometryIterator_;
  typedef GeometryIterator<OBS>       ObservationIterator_;
  typedef IncrementSet<MODEL>         IncrementSet_;
  typedef ObsErrors<OBS>              ObsErrors_;
  typedef ObsSpace<OBS>               ObsSpace_;
  typedef ObsSpaces<OBS>              ObsSpaces_;
  typedef ObsVector<OBS>              ObsVector_;
  typedef StateSet<MODEL>             StateSet_;

 public:
  static const std::string classname() {return "oops::SequentialEnsembleSolver";}

  SequentialEnsembleSolver(ObsSpaces_ & obspaces, const Geometry_ & geometry,
                      const eckit::Configuration & config, size_t nens, const StateSet_ & xbmean,
                      const Variables & incvars);

  virtual ~SequentialEnsembleSolver() = default;

  /// update background ensemble \p bg to analysis ensemble \p an for all points on this PE
  void measurementUpdate(const IncrementSet_ & bg,
                         IncrementSet_ & an) override;

  // IGNORE THIS, hack for compiling purposes only
  void measurementUpdate(const Eigen::VectorXd &,
                         const ObsErrors_ &,
                         const Departures_ &,
                         const IncrementSet_ &,
                         const GeometryIterator_ &,
                         IncrementSet_ &)  override {
    ABORT("SequentialEnsembleSolver::measurementUpdate(6) should never be called");
  }

 protected:
  /// Update the observation prior ensemble for a single observation, this is the
  /// core of the sequential DA update and implementation varies based on
  /// the flavor of sequential DA being used.
  ///
  /// \param yb_k Observation ensemble perturbations priors
  /// \param omb_k Observation minus background mean (innovation) for the k-th observation
  /// \param oberr_variance_k Observation error variance for the k-th observation
  /// \param delta_y_k Output observation ensemble increment for the k-th observation
  ///                  to be updated in-place by the derived class implementation
virtual void obsEnsembleUpdate(const Eigen::VectorXf & yb_k,
                               const double omb_k,
                               const double oberr_variance_k,
                               Eigen::VectorXf & delta_y_k) = 0;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
SequentialEnsembleSolver<MODEL, OBS>::SequentialEnsembleSolver(ObsSpaces_ & obspaces,
                                                     const Geometry_ & geometry,
                                                     const eckit::Configuration & config,
                                                     size_t nens,
                                                     const StateSet_ & xbmean,
                                                     const Variables & incvars)
      : LocalEnsembleSolver<MODEL, OBS>(obspaces, geometry, config, nens, xbmean, incvars)
  { }

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void SequentialEnsembleSolver<MODEL, OBS>::measurementUpdate(const IncrementSet_ & bkg_pert,
                                                             IncrementSet_ & ana_pert) {
  // NOTE: assumes the obspaces share this communicator; oops::ObsSpace has no
  // portable spatial-comm accessor to assert it.
  const eckit::mpi::Comm & comm = this->geometry_.getComm();

  // initialize analysis ensemble to background ensemble
  ana_pert = bkg_pert;

  // prep the local variables that are needed
  Departures_ mask(this->obspaces_);
  mask.ones();
  this->applyAssimilatedMask(mask);
  Eigen::VectorXd omb = this->omb_.packEigen(mask);
  Eigen::MatrixXf yb = (this->Yb_)->packEigen(mask);
  Eigen::VectorXd invVarR = this->invVarR_->packEigen(mask);
  std::vector<double> maskVec;  // 1.0 where valid, some other crazy number where masked out
  for (size_t i = 0; i < mask.size(); ++i) {
    std::vector<double> maskVecSub;
    mask[i].serialize(maskVecSub);
    maskVec.insert(maskVec.end(), maskVecSub.begin(), maskVecSub.end());
  }

  // sanity check, we have more than 1 ensemble member, right!?
  if (yb.rows() <= 1) {
    ABORT("SequentialEnsembleSolver::measurementUpdate requires more than 1 ensemble member");
  }

  // create a copy of the locations for the local obs.
  // I COULD have just used the iterator in the loops later, but the iterator
  // gets complicated when there are multiple obs spaces and obs spaces with
  // more than one variable, so i'm just throwing all that annoyance together
  // here instead.
  std::vector<eckit::geometry::Point3> locations;
  locations.reserve(omb.size());
  std::vector<double>::iterator maskVecItr = maskVec.begin();
  for (size_t obSpaceIdx = 0; obSpaceIdx < this->obspaces_.size(); obSpaceIdx++) {
    const ObsSpace_ & obspace = this->obspaces_[obSpaceIdx];
    // note, max(1,...) is used because some toy models return 0 variables (e.g. L95)
    const size_t nObspaceVars = std::max(static_cast<size_t>(1),
                                         obspace.obsvariables().size());
    for (ObservationIterator_ iter = obspace.begin(); iter != obspace.end(); ++iter) {
      if (*maskVecItr == 1.0) {
        for (size_t jj = 0; jj < nObspaceVars; jj++) {
          locations.push_back(*iter);
        }
      }
      ++maskVecItr;
    }
  }
  ASSERT(static_cast<Eigen::Index>(locations.size()) == omb.size());

  // create a distributed list of which PE owns which observations
  const int nobs_local = omb.size();
  int nobs_global;
  comm.allReduce(nobs_local, nobs_global, eckit::mpi::sum());
  std::vector<int> obOwnerRank(nobs_global);
  {
    // TODO(Travis) change this to a round-robin distribution (or get the indexes
    // from ioda?) Right now this processes the obs in PE order, which is not ideal as
    // it will lead to changes in answers with changes in PE count (since
    // sequential EnKF is sensitive to order the obs are assimilated in).
    std::vector<int> obOwnerRank_local(nobs_local, comm.rank());
    std::vector<int> peRecvCount(comm.size());
    std::vector<int> peDispls(comm.size());
    comm.allGather(nobs_local, peRecvCount.begin(), peRecvCount.end());
    std::exclusive_scan(peRecvCount.begin(), peRecvCount.end(), peDispls.begin(), 0);
    comm.allGatherv(obOwnerRank_local.begin(), obOwnerRank_local.end(), obOwnerRank.begin(),
                    peRecvCount.data(), peDispls.data());
  }

  // iterate globally over each observation sequentially
  Eigen::Index nextLocalObIdx = 0;  // next ob on THIS PE to assimilate when it's our turn
  for (Eigen::Index kk = 0; kk < nobs_global; ++kk)  {
    // determine which PE this should be coming from
    size_t obSourcePE = obOwnerRank[kk];
    bool obOwnedByThisPE = obSourcePE == comm.rank();

    // observation prior ensemble and increments, for the current observation
    // (to be calculated by the owner PE and then broadcast to all PEs)
    Eigen::VectorXf yb_k(yb.rows());
    Eigen::VectorXf delta_y_k(yb_k.size());
    eckit::geometry::Point3 y_location;

    // calculate the observation increment (if PE owns this ob)
    if (obOwnedByThisPE) {
      if (nextLocalObIdx >= omb.size()) {
        // we shouldn't get here, handle this
        ABORT("Error on PE running past local obs array");
      }

      // get all the needed info for this obs (if PE owns this ob)
      // note, sanity check is done to ob error variance to make sure it is a
      // non-zero positive number, just in case.
      double omb_k = omb(nextLocalObIdx);
      double oberr_k = std::max(1.0 / invVarR(nextLocalObIdx), 1e-12);
      yb_k = yb.col(nextLocalObIdx);
      y_location = locations[nextLocalObIdx];

      // calculate the observation ensemble increments (if PE owns this ob)
      this->obsEnsembleUpdate(yb_k, omb_k, oberr_k, delta_y_k);
    }

    // broadcast the packed data to all PEs
    // TODO(Travis) check if any communication time can be saved by packing this
    // all into a single vector first
    comm.broadcast(yb_k.data(), yb_k.size(), obSourcePE);
    comm.broadcast(delta_y_k.data(), delta_y_k.size(), obSourcePE);
    comm.broadcast(y_location.data(), 3, obSourcePE);

    // compute deviations for observation ensemble, and squared norm, for use later.
    // Explicit recentering is required: yb_k is not guaranteed to have zero mean at
    // this point (verified empirically — omitting the subtraction changes analysis
    // answers in the L95 and QG tests).
    Eigen::VectorXf yb_k_dev = yb_k.array() - yb_k.mean();
    double yb_k_dev2 = yb_k_dev.squaredNorm();

    // if there is no spread in the obs ensemble, skip it. otherwise, update the
    // subsequent obs ensembles and state ensemble
    if (yb_k_dev2 > 0.0) {
      // iterate over all unused observation prior ensembles on this PE, to update them.
      // skip the current observation (jj == nextLocalObIdx) on the owning PE since it
      // has already been assimilated and updating it would be wasteful.
      const Eigen::Index jj_start = nextLocalObIdx + (obOwnedByThisPE ? 1 : 0);
      for (Eigen::Index jj = jj_start; jj < nobs_local; ++jj) {
        // calculate localization between the two obs
        double localization = this->obsloc().computeLocalization(locations[jj], y_location);

        // if localization > 0.0, update the observation prior ensemble
        if (localization > 0.0)  {
          Eigen::VectorXf yb_j = yb.col(jj);
          Eigen::VectorXf delta_y_j(delta_y_k.size());

          // calculate regression between the two observation ensembles
          double beta = yb_k_dev.dot(yb_j) / yb_k_dev2;

          // update the observation prior ensembles, and omb.
          // ensure mean of the delta_y_j is zero, since omb is updated separately with the mean
          delta_y_j = localization * beta * delta_y_k;
          const double delta_mean = delta_y_j.mean();
          delta_y_j.array() -= delta_mean;
          yb_j += delta_y_j;
          yb.col(jj) = yb_j;
          omb(jj) -= delta_mean;
        }
      }

      // state ensemble update
      for (GeometryIterator_ geomItr = this->geometry_.begin();
          geomItr != this->geometry_.end(); ++geomItr)  {
        // calculate localization between this grid point and the observation
        // location, and skip update if localization <= 0.0
        double localization = this->obsloc().computeLocalization(*geomItr, y_location);
        if (localization <= 0.0) continue;

        // for each time
        for (size_t itime = 0; itime < ana_pert.time_size(); ++itime) {
          // get the local state ensemble at this grid point
          // TODO(someone) this pack/set is inefficient, and (at least for QG) is
          // a main source of unnecessary overhead.
          Eigen::MatrixXd Xa;
          ana_pert.packEigen(Xa, geomItr, itime);

          // calculate the regression between observation and state ensembles.
          Eigen::VectorXd beta = (Xa * yb_k_dev.cast<double>()) / yb_k_dev2;

          // update the local state ensemble at this grid point
          Xa += localization * beta * delta_y_k.transpose().cast<double>();
          ana_pert.setEigen(Xa, geomItr, itime);
        }
      }  // end of state ensemble update
    }

    // increment the observation index
    if (obOwnedByThisPE) {
      ++nextLocalObIdx;
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_ASSIMILATION_SEQUENTIALENSEMBLESOLVER_H_
