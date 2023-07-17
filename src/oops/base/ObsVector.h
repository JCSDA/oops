/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSVECTOR_H_
#define OOPS_BASE_OBSVECTOR_H_

#include <math.h>
#include <memory>
#include <ostream>
#include <string>
#include <utility>

#include "oops/interface/ObsVector.h"
#include "oops/util/gatherPrint.h"

namespace oops {

template<typename OBS> class ObsSpace;

// -----------------------------------------------------------------------------
/// \brief ObsVector class used in oops; subclass of interface class interface::ObsVector.
///
/// \details
/// Handles additional MPI communicator parameter \p commTime_ in the constructors
/// (for MPI distribution in time, used in oops for 4DEnVar and weak-constraint 4DVar).
/// Adds communication through time to norm and print methods.
// -----------------------------------------------------------------------------

template <typename OBS>
class ObsVector : public interface::ObsVector<OBS>  {
  typedef typename OBS::ObsVector             ObsVector_;

 public:
  static const std::string classname() {return "oops::ObsVector";}

  /// Creates vector from the simulated variable in the \p obsspace. If \p name is specified,
  /// reads the specified \p name variable from \p obsspace. Otherwise, zero vector is created.
  explicit ObsVector(const ObsSpace<OBS> & obsspace, const std::string name = "");

  /// Wraps an existing ObsVector_.
  ///
  /// \param obsvector  The vector to wrap.
  /// \param timeComm   Time communicator.
  ObsVector(std::unique_ptr<ObsVector_> obsvector, const eckit::mpi::Comm &timeComm);
  /// Copy constructor
  ObsVector(const ObsVector &);

  /// Use assignment operator (const ObsDataVector<OBS, float> &) from the base class
  using interface::ObsVector<OBS>::operator=;

  /// Return the dot product between this ObsVector and \p other ObsVector
  double dot_product_with(const ObsVector & other) const;
  /// Return this ObsVector rms
  double rms() const;
  /// Number of non-masked out observations (across all MPI tasks)
  unsigned int nobs() const;

 private:
  void print(std::ostream &) const;
  const eckit::mpi::Comm * commTime_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsVector<OBS>::ObsVector(const ObsSpace<OBS> & os, const std::string name)
  : interface::ObsVector<OBS>(os, name), commTime_(&os.timeComm()) {}

// -----------------------------------------------------------------------------

template <typename OBS>
ObsVector<OBS>::ObsVector(std::unique_ptr<ObsVector_> obsvector, const eckit::mpi::Comm &timeComm)
  : interface::ObsVector<OBS>(std::move(obsvector)), commTime_(&timeComm) {}

// -----------------------------------------------------------------------------

template <typename OBS>
ObsVector<OBS>::ObsVector(const ObsVector & other): interface::ObsVector<OBS>(other),
  commTime_(other.commTime_) {}

// -----------------------------------------------------------------------------

template <typename OBS>
double ObsVector<OBS>::dot_product_with(const ObsVector & other) const {
  double zz = interface::ObsVector<OBS>::dot_product_with(other);
  commTime_->allReduceInPlace(zz, eckit::mpi::Operation::SUM);
  return zz;
}

// -----------------------------------------------------------------------------

template <typename OBS>
double ObsVector<OBS>::rms() const {
  double zz = interface::ObsVector<OBS>::rms();
  size_t ntot = interface::ObsVector<OBS>::nobs();
  double zzz = zz * zz * static_cast<double>(ntot);
  commTime_->allReduceInPlace(zzz, eckit::mpi::Operation::SUM);
  ntot = nobs();
  if (ntot > 0) {
    zzz /= static_cast<double>(ntot);
    zz = std::sqrt(zzz);
  } else {
    zz = 0.0;
  }
  return zz;
}

// -----------------------------------------------------------------------------

template <typename OBS>
unsigned int ObsVector<OBS>::nobs() const {
  int nobs = interface::ObsVector<OBS>::nobs();
  commTime_->allReduceInPlace(nobs, eckit::mpi::Operation::SUM);
  return nobs;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsVector<OBS>::print(std::ostream & os) const {
  if (commTime_->size() > 1) {
    util::gatherPrint(os, this->obsvector(), *commTime_);
  } else {
    os << this->obsvector();
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSVECTOR_H_
