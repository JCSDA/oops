/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_LOCALIZATION_H_
#define OOPS_BASE_LOCALIZATION_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/generic/LocalizationBase.h"
#include "oops/mpi/mpi.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

/// \brief Abstract model-space localization class used by high level algorithms
/// and applications.
///
/// Note: to see methods that need to be implemented in a generic Localization
/// implementation, see LocalizationBase class in generic/LocalizationBase.h.
/// To see methods that need to be implemented in a MODEL-specific Localization
/// implementation, see interface::LocalizationBase class in
/// interface/LocalizationBase.h.
template <typename MODEL>
class Localization : public util::Printable,
                     private boost::noncopyable,
                     private util::ObjectCounter<Localization<MODEL> > {
  typedef Geometry<MODEL>   Geometry_;
  typedef Increment<MODEL>  Increment_;
  typedef LocalizationBase<MODEL>   LocBase_;

 public:
  static const std::string classname() {return "oops::Localization";}

  /// Set up Localization for \p geometry, configured with \p conf
  Localization(const Geometry_ & geometry,
               const oops::Variables & incVars,
               const eckit::Configuration & conf);
  ~Localization();

  /// Randomize \p dx and apply 4D localization. All 3D blocks of the 4D localization
  /// matrix are the same (and defined by 3D localization loc_)
  virtual void randomize(Increment_ & dx) const;
  /// Apply 4D localization. All 3D blocks of the 4D localization matrix are the same
  /// (and defined by 3D localization loc_)
  virtual void multiply(Increment_ & dx) const;

 private:
  /// Print, used in logging
  void print(std::ostream &) const override;

  std::unique_ptr<util::Timer> timeConstr_;
  /// Pointer to the Localization implementation
  std::unique_ptr<LocBase_> loc_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
Localization<MODEL>::Localization(const Geometry_ & geometry,
                                  const oops::Variables & incVars,
                                  const eckit::Configuration & conf)
  : timeConstr_(new util::Timer(classname(), "Localization")),
    loc_(LocalizationFactory<MODEL>::create(geometry, incVars, conf))
{
  Log::trace() << "Localization<MODEL>::Localization done" << std::endl;
  timeConstr_.reset();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Localization<MODEL>::~Localization() {
  Log::trace() << "Localization<MODEL>::~Localization starting" << std::endl;
  util::Timer timer(classname(), "~Localization");
  loc_.reset();
  Log::trace() << "Localization<MODEL>::~Localization done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void Localization<MODEL>::randomize(Increment_ & dx) const {
  Log::trace() << "Localization<MODEL>::randomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");
  const eckit::mpi::Comm & comm = dx.timeComm();
  static int tag = 23456;
  size_t nslots = comm.size();
  int mytime = comm.rank();

  if (mytime > 0) {
    util::DateTime dt = dx.validTime();   // Save original time value
    dx.zero();
    oops::mpi::receive(comm, dx, 0, tag);
    dx.updateTime(dt - dx.validTime());  // Set time back to original value
  } else {
    // Apply 3D localization
    loc_->randomize(dx);

    // Copy result to all timeslots
    for (size_t jj = 1; jj < nslots; ++jj) {
      oops::mpi::send(comm, dx, jj, tag);
    }
  }
  ++tag;
  Log::trace() << "Localization<MODEL>::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void Localization<MODEL>::multiply(Increment_ & dx) const {
  Log::trace() << "Localization<MODEL>::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  const eckit::mpi::Comm & comm = dx.timeComm();
  static int tag = 23456;
  size_t nslots = comm.size();
  int mytime = comm.rank();

  // Use Mark Buehner's trick to save CPU when applying the same 3D localization for all
  // 3D blocks of the 4D localization matrix:
  // L_4D = ( L_3D L_3D L_3D ) = ( Id ) L_3D ( Id Id Id )
  //        ( L_3D L_3D L_3D )   ( Id )
  //        ( L_3D L_3D L_3D )   ( Id )
  // so if :
  // x_4D = ( x_1 )
  //        ( x_2 )
  //        ( x_3 )
  // then:
  // L_4D x_4D = (Id) L_3D (Id Id Id) (x_1) = (Id) L_3D (x_1+x_2+x_3) = (L_3D ( x_1 + x_2 + x_3 ))
  //             (Id)                 (x_2)   (Id)                      (L_3D ( x_1 + x_2 + x_3 ))
  //             (Id)                 (x_3)   (Id)                      (L_3D ( x_1 + x_2 + x_3 ))
  // Reference in section 3.4.2. of https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.2325.
  if (mytime > 0) {
    util::DateTime dt = dx.validTime();   // Save original time value
    oops::mpi::send(comm, dx, 0, tag);
    dx.zero();
    oops::mpi::receive(comm, dx, 0, tag);
    dx.updateTime(dt - dx.validTime());  // Set time back to original value
  } else {
    // Sum over timeslots
    for (size_t jj = 1; jj < nslots; ++jj) {
      Increment_ dxtmp(dx);
      oops::mpi::receive(comm, dxtmp, jj, tag);
      dx.axpy(1.0, dxtmp, false);
    }

    // Apply 3D localization
    loc_->multiply(dx);

    // Copy result to all timeslots
    for (size_t jj = 1; jj < nslots; ++jj) {
      oops::mpi::send(comm, dx, jj, tag);
    }
  }
  ++tag;

  Log::trace() << "Localization<MODEL>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void Localization<MODEL>::print(std::ostream & os) const {
  Log::trace() << "Localization<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *loc_;
  Log::trace() << "Localization<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_LOCALIZATION_H_
