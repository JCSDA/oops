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

#ifndef OOPS_INTERFACE_LOCATIONS_H_
#define OOPS_INTERFACE_LOCATIONS_H_

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
/// \brief Locations of observations for observation operator
template <typename OBS>
class Locations : public util::Printable,
                  private util::ObjectCounter<Locations<OBS> > {
  typedef typename OBS::Locations             Locations_;
 public:
  static const std::string classname() {return "oops::Locations";}

  /// Constructor from the pointer returned by ObsOperator::locations()
  explicit Locations(std::unique_ptr<Locations_>);
  /// Constructor used in tests
  Locations(const eckit::Configuration &, const eckit::mpi::Comm &);

  /// Destructor and copy/move constructor and assignments
  ~Locations();
  Locations(const Locations &) = delete;
  Locations(Locations &&);
  Locations & operator=(const Locations &) = delete;
  Locations & operator=(Locations &&);

  const std::vector<double> & latitudes() const;
  const std::vector<double> & longitudes() const;
  const std::vector<util::DateTime> & times() const;

  /// Interfacing
  const Locations_ & locations() const {return *locs_;}

 private:
  void print(std::ostream &) const;
  std::unique_ptr<const Locations_> locs_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
Locations<OBS>::Locations(std::unique_ptr<Locations_> locs) : locs_(std::move(locs)) {
  Log::trace() << "Locations<OBS>::Locations constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
Locations<OBS>::Locations(const eckit::Configuration & conf, const eckit::mpi::Comm & comm) {
  Log::trace() << "Locations<OBS>::Locations starting" << std::endl;
  util::Timer timer(classname(), "Locations");
  locs_.reset(new Locations_(conf, comm));
  Log::trace() << "Locations<OBS>::Locations done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
Locations<OBS>::~Locations() {
  Log::trace() << "Locations<OBS>::~Locations starting" << std::endl;
  util::Timer timer(classname(), "~Locations");
  locs_.reset();
  Log::trace() << "Locations<OBS>::~Locations done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
Locations<OBS>::Locations(Locations && other): locs_(std::move(other.locs_)) {
  util::Timer timer(classname(), "Locations");
  Log::trace() << "Locations<OBS> moved" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
Locations<OBS> & Locations<OBS>::operator=(Locations<OBS> && other) {
  Log::trace() << "Locations<OBS>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  locs_ = std::move(other.locs_);
  Log::trace() << "Locations<OBS>::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template <typename OBS>
const std::vector<double> & Locations<OBS>::latitudes() const {
  Log::trace() << "Locations<OBS>::latitudes starting" << std::endl;
  util::Timer timer(classname(), "latitudes");
  return locs_->latitudes();
}

// -----------------------------------------------------------------------------

template <typename OBS>
const std::vector<double> & Locations<OBS>::longitudes() const {
  Log::trace() << "Locations<OBS>::longitudes starting" << std::endl;
  util::Timer timer(classname(), "longitudes");
  return locs_->longitudes();
}

// -----------------------------------------------------------------------------

template <typename OBS>
const std::vector<util::DateTime> & Locations<OBS>::times() const {
  Log::trace() << "Locations<OBS>::times starting" << std::endl;
  util::Timer timer(classname(), "times");
  return locs_->times();
}

// -----------------------------------------------------------------------------

template<typename OBS>
void Locations<OBS>::print(std::ostream & os) const {
  Log::trace() << "Locations<OBS>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *locs_;
  Log::trace() << "Locations<OBS>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_LOCATIONS_H_
