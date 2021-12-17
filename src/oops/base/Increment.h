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

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/interface/Increment.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/gatherPrint.h"
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

 public:
  /// Constructor for specified \p geometry, with \p variables, valid on \p date
  Increment(const Geometry_ & geometry, const Variables & variables, const util::DateTime & date);
  /// Copies \p other increment, changing its resolution to \p geometry
  Increment(const Geometry_ & geometry, const Increment & other);
  /// Creates Increment with the same geometry and variables as \p other.
  /// Copies \p other if \p copy is true, otherwise creates zero increment
  Increment(const Increment & other, const bool copy = true);

  /// Accessor to the time communicator
  const eckit::mpi::Comm & timeComm() const {return *timeComm_;}

  /// Shift forward in time by \p dt
  void shift_forward(const util::DateTime & dt);
  /// Shift backward in time by \p dt
  void shift_backward(const util::DateTime & dt);

  /// Set ATLAS fieldset associated with this Increment internally
  void toAtlas();
  /// Allow to access base class's method as well
  using interface::Increment<MODEL>::toAtlas;
  /// Accessors to the ATLAS fieldset
  atlas::FieldSet & atlas() {return atlasFieldSet_;}
  const atlas::FieldSet & atlas() const {return atlasFieldSet_;}

  /// dot product with the \p other increment
  double dot_product_with(const Increment & other) const;
  /// Norm for diagnostics
  double norm() const;

 private:
  void print(std::ostream &) const override;

  const eckit::mpi::Comm * timeComm_;  /// pointer to the MPI communicator in time
  atlas::FieldSet atlasFieldSet_;      /// Atlas fields associated with this Increment
};

// -----------------------------------------------------------------------------

template <typename MODEL>
Increment<MODEL>::Increment(const Geometry_ & geometry, const Variables & variables,
                            const util::DateTime & date):
  interface::Increment<MODEL>(geometry, variables, date),
  timeComm_(&geometry.timeComm())
{}

// -----------------------------------------------------------------------------

template <typename MODEL>
Increment<MODEL>::Increment(const Geometry_ & geometry, const Increment & other):
  interface::Increment<MODEL>(geometry, other),
  timeComm_(other.timeComm_)
{}

// -----------------------------------------------------------------------------

template <typename MODEL>
Increment<MODEL>::Increment(const Increment & other, const bool copy):
  interface::Increment<MODEL>(other, copy),
  timeComm_(other.timeComm_)
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
double Increment<MODEL>::dot_product_with(const Increment & dx) const {
  double zz = interface::Increment<MODEL>::dot_product_with(dx);
  timeComm_->allReduceInPlace(zz, eckit::mpi::Operation::SUM);
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double Increment<MODEL>::norm() const {
  double zz = interface::Increment<MODEL>::norm();
  zz *= zz;
  timeComm_->allReduceInPlace(zz, eckit::mpi::Operation::SUM);
  zz = sqrt(zz);
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::shift_forward(const util::DateTime & begin) {
  Log::trace() << "Increment<MODEL>::Increment shift_forward starting" << std::endl;
  static int tag = 159357;
  size_t mytime = timeComm_->rank();

// Send values of M.dx_i at end of my subwindow to next subwindow
  if (mytime + 1 < timeComm_->size()) {
    oops::mpi::send(*timeComm_, *this, mytime+1, tag);
  }

// Receive values at beginning of my subwindow from previous subwindow
  if (mytime > 0) {
    oops::mpi::receive(*timeComm_, *this, mytime-1, tag);
  } else {
    this->zero(begin);
  }

  ++tag;
  Log::trace() << "Increment<MODEL>::Increment shift_forward done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::shift_backward(const util::DateTime & end) {
  Log::trace() << "Increment<MODEL>::Increment shift_backward starting" << std::endl;
  static int tag = 30951;
  size_t mytime = timeComm_->rank();

// Send values of dx_i at start of my subwindow to previous subwindow
  if (mytime > 0) {
    oops::mpi::send(*timeComm_, *this, mytime-1, tag);
  }

// Receive values at end of my subwindow from next subwindow
  if (mytime + 1 < timeComm_->size()) {
    oops::mpi::receive(*timeComm_, *this, mytime+1, tag);
  } else {
    this->zero(end);
  }

  ++tag;
  Log::trace() << "Increment<MODEL>::Increment shift_backward done" << std::endl;
}


// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment<MODEL>::toAtlas() {
  interface::Increment<MODEL>::setAtlas(&atlasFieldSet_);
  interface::Increment<MODEL>::toAtlas(&atlasFieldSet_);
  this->increment_.reset();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::print(std::ostream & os) const {
  if (timeComm_->size() > 1) {
    gatherPrint(os, this->increment(), *timeComm_);
  } else {
    os << this->increment();
  }
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


}  // namespace oops

#endif  // OOPS_BASE_INCREMENT_H_
