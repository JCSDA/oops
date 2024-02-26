/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_SAMPLEDLOCATIONS_H_
#define OOPS_INTERFACE_SAMPLEDLOCATIONS_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// \brief A set of interpolation paths sampling the observation locations.
///
/// The location of a single observation may be an extended region of space and be sampled (crossed)
/// by more than one path. Conversely, some paths may sample locations of multiple observations.
///
/// \note Currently each interpolation path is assumed to be vertical. There are plans to add
/// support for slanted paths in the future.
template <typename OBS>
class SampledLocations : public util::Printable,
                         private util::ObjectCounter<SampledLocations<OBS>> {
  typedef typename OBS::SampledLocations SampledLocations_;
 public:
  static const std::string classname() { return "oops::SampledLocations"; }

  /// Wrap a pointer to an implementation of this interface
  explicit SampledLocations(std::unique_ptr<SampledLocations_>);
  /// Constructor used in tests
  SampledLocations(const eckit::Configuration &, const eckit::mpi::Comm &);

  /// Destructor and copy/move constructor and assignments
  ~SampledLocations();
  SampledLocations(const SampledLocations &) = delete;
  SampledLocations(SampledLocations &&);
  SampledLocations & operator=(const SampledLocations &) = delete;
  SampledLocations & operator=(SampledLocations &&);

  /// Latitudes of the (vertical) interpolation paths.
  const std::vector<double> & latitudes() const;
  /// Longitudes of the (vertical) interpolation paths.
  const std::vector<double> & longitudes() const;
  /// Interpolation times.
  const std::vector<util::DateTime> & times() const;

  /// Interfacing
  const SampledLocations_ & sampledLocations() const {return *sampledLocations_;}

 private:
  void print(std::ostream &) const;
  std::unique_ptr<const SampledLocations_> sampledLocations_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
SampledLocations<OBS>::SampledLocations(std::unique_ptr<SampledLocations_> sampledLocations)
  : sampledLocations_(std::move(sampledLocations)) {
  Log::trace() << "SampledLocations<OBS>::SampledLocations constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
SampledLocations<OBS>::SampledLocations(const eckit::Configuration & conf,
                                        const eckit::mpi::Comm & comm) {
  Log::trace() << "SampledLocations<OBS>::SampledLocations starting" << std::endl;
  util::Timer timer(classname(), "SampledLocations");
  sampledLocations_.reset(new SampledLocations_(conf, comm));
  Log::trace() << "SampledLocations<OBS>::SampledLocations done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
SampledLocations<OBS>::~SampledLocations() {
  Log::trace() << "SampledLocations<OBS>::~SampledLocations starting" << std::endl;
  util::Timer timer(classname(), "~SampledLocations");
  sampledLocations_.reset();
  Log::trace() << "SampledLocations<OBS>::~SampledLocations done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
SampledLocations<OBS>::SampledLocations(SampledLocations && other):
    sampledLocations_(std::move(other.sampledLocations_)) {
  util::Timer timer(classname(), "SampledLocations");
  Log::trace() << "SampledLocations<OBS> moved" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
SampledLocations<OBS> & SampledLocations<OBS>::operator=(SampledLocations<OBS> && other) {
  Log::trace() << "SampledLocations<OBS>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  sampledLocations_ = std::move(other.sampledLocations_);
  Log::trace() << "SampledLocations<OBS>::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template <typename OBS>
const std::vector<double> & SampledLocations<OBS>::latitudes() const {
  Log::trace() << "SampledLocations<OBS>::latitudes starting" << std::endl;
  util::Timer timer(classname(), "latitudes");
  return sampledLocations_->latitudes();
}

// -----------------------------------------------------------------------------

template <typename OBS>
const std::vector<double> & SampledLocations<OBS>::longitudes() const {
  Log::trace() << "SampledLocations<OBS>::longitudes starting" << std::endl;
  util::Timer timer(classname(), "longitudes");
  return sampledLocations_->longitudes();
}

// -----------------------------------------------------------------------------

template <typename OBS>
const std::vector<util::DateTime> & SampledLocations<OBS>::times() const {
  Log::trace() << "SampledLocations<OBS>::times starting" << std::endl;
  util::Timer timer(classname(), "times");
  return sampledLocations_->times();
}

// -----------------------------------------------------------------------------

template<typename OBS>
void SampledLocations<OBS>::print(std::ostream & os) const {
  Log::trace() << "SampledLocations<OBS>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *sampledLocations_;
  Log::trace() << "SampledLocations<OBS>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_SAMPLEDLOCATIONS_H_
