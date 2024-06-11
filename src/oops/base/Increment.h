/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_INCREMENT_H_
#define OOPS_BASE_INCREMENT_H_

#include <memory>
#include<string>
#include<vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/interface/Increment.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Timer.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Increment class used in oops
///
/// \details
/// Adds extra methods that do not need to be implemented in the model implementations:
/// - timeComm()  (accessor to the MPI communicator in time - collection of processes
///                holding the data needed to represent the state in a particular region
///                of space X_i and throughout the whole time interval for which DA is done)
/// - shift_forward
/// - shift_backward
/// - toAtlas, atlas
///
/// Adds communication through time to the following Increment methods:
/// - dot_product_with
/// - norm
/// - print

template <typename MODEL>
class Increment : public interface::Increment<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;

 public:
  /// Constructor for specified \p geometry, with \p variables, valid on \p date
  Increment(const Geometry_ & geometry, const Variables & variables, const util::DateTime & date);
  /// Copies \p other increment, changing its resolution to \p geometry
  /// Uses adjoint resolution change if \p ad is true, though still to the target \p geometry
  Increment(const Geometry_ & geometry, const Increment & other, const bool ad = false);
  /// Creates Increment with the same geometry and variables as \p other.
  /// Copies \p other if \p copy is true, otherwise creates zero increment
  Increment(const Increment & other, const bool copy = true);

  Increment & operator=(const Increment &);

  /// Accessor to geometry associated with this Increment
  const Geometry_ & geometry() const {return resol_;}

  void transfer_from_state(const State_ &);

  /// Accessors to the FieldSet3D
  const FieldSet3D & fieldSet() const;
  FieldSet3D & fieldSet();
  void synchronizeFields();

  /// Compute root-mean-square by variable by level
  /// For preconditioning HybridLinearModel coefficient calculation
  std::vector<double> rmsByVariableByLevel(const Variable &, const bool) const;

 private:
  const Geometry_ & resol_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
Increment<MODEL>::Increment(const Geometry_ & geometry, const Variables & variables,
                            const util::DateTime & date):
  interface::Increment<MODEL>(geometry, variables, date), resol_(geometry)
{}

// -----------------------------------------------------------------------------

template <typename MODEL>
Increment<MODEL>::Increment(const Geometry_ & geometry, const Increment & other, const bool ad):
  interface::Increment<MODEL>(geometry, other, ad), resol_(geometry)
{}

// -----------------------------------------------------------------------------

template <typename MODEL>
Increment<MODEL>::Increment(const Increment & other, const bool copy):
  interface::Increment<MODEL>(other, copy), resol_(other.resol_)
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> & Increment<MODEL>::operator=(const Increment & rhs) {
  ASSERT(resol_ == rhs.resol_);
  interface::Increment<MODEL>::operator=(rhs);
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
const FieldSet3D & Increment<MODEL>::fieldSet() const {
  if (!interface::Increment<MODEL>::fset_) {
    interface::Increment<MODEL>::fset_.reset(new FieldSet3D(this->validTime(), resol_.getComm()));
  }
  if (interface::Increment<MODEL>::fset_->empty()) {
    this->toFieldSet(interface::Increment<MODEL>::fset_->fieldSet());
    for (const auto & field : *interface::Increment<MODEL>::fset_) {
      ASSERT_MSG(field.rank() == 2,
                 "OOPS expects the model's Increment::toFieldSet method to return rank-2 fields,"
                 " but field " + field.name() + " has rank = " + std::to_string(field.rank()));
    }
  }
  return *interface::Increment<MODEL>::fset_;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
FieldSet3D & Increment<MODEL>::fieldSet() {
  if (!interface::Increment<MODEL>::fset_) {
    interface::Increment<MODEL>::fset_.reset(new FieldSet3D(this->validTime(), resol_.getComm()));
  }
  if (interface::Increment<MODEL>::fset_->empty()) {
    this->toFieldSet(interface::Increment<MODEL>::fset_->fieldSet());
    for (const auto & field : *interface::Increment<MODEL>::fset_) {
      ASSERT_MSG(field.rank() == 2,
                 "OOPS expects the model's Increment::toFieldSet method to return rank-2 fields,"
                 " but field " + field.name() + " has rank = " + std::to_string(field.rank()));
    }
  }
  return *interface::Increment<MODEL>::fset_;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::synchronizeFields() {
  // TODO(JEDI core team): remove this method when accessors are fully implemented
  if (interface::Increment<MODEL>::fset_) {
    ASSERT(!interface::Increment<MODEL>::fset_->empty());
    this->fromFieldSet(interface::Increment<MODEL>::fset_->fieldSet());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::transfer_from_state(const State_ & xx) {
  // TODO(JEDI core team): this is inneficient, we should recode this method using fieldSets
  // when those are fully implemented
  State_ zz(xx);
  zz.zero();
  this->diff(xx, zz);
}

// -----------------------------------------------------------------------------

/// Add on \p dx incrment to model state \p xx
template <typename MODEL>
State<MODEL> & operator+=(State<MODEL> & xx, const Increment<MODEL> & dx) {
  Log::trace() << "operator+=(State, Increment) starting" << std::endl;
  util::Timer timer("oops::Increment", "operator+=(State, Increment)");
  xx.state() += dx.increment();
  Log::trace() << "operator+=(State, Increment) done" << std::endl;
  return xx;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
std::vector<double> Increment<MODEL>::rmsByVariableByLevel(const Variable & var,
                                                           const bool global) const {
  Log::trace() << "Increment<MODEL>::rmsByVariableByLevel starting" << std::endl;
  util::Timer timer("oops::Increment", "rmsByVariableByLevel");
  if (resol_.fields().empty()) {
    ABORT("Increment<MODEL>::rmsByVariableByLevel: requires Atlas interface");
  }
  const auto ownedView = atlas::array::make_view<int, 2>(resol_.fields()["owned"]);
  const auto fieldView = atlas::array::make_view<double, 2>(fieldSet()[var.name()]);
  std::vector<double> rms(fieldView.shape(1), 0.0);
  for (atlas::idx_t k = 0; k < fieldView.shape(1); ++k) {
    size_t nOwned = 0;
    for (atlas::idx_t i = 0; i < fieldView.shape(0); ++i) {
      if (ownedView(i, 0) > 0) {
        ++nOwned;
        rms[k] += fieldView(i, k) * fieldView(i, k);
      }
    }
    if (global) {
      resol_.getComm().allReduceInPlace(nOwned, eckit::mpi::sum());
      resol_.getComm().allReduceInPlace(rms[k], eckit::mpi::sum());
    }
    rms[k] = sqrt(rms[k] / nOwned);
  }
  return rms;
  Log::trace() << "Increment<MODEL>::rmsByVariableByLevel done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_INCREMENT_H_
