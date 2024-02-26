/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2021-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_LOCSL95_H_
#define LORENZ95_LOCSL95_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace lorenz95 {

/// Implements the oops SampledLocations interface for the L95 obs type.
class LocsL95 : public util::Printable,
                private util::ObjectCounter<LocsL95> {
 public:
  static const std::string classname() {return "lorenz95::LocsL95";}

  LocsL95(const std::vector<double> &, const std::vector<util::DateTime> &);
  LocsL95(const eckit::Configuration &, const eckit::mpi::Comm &);

  size_t size() const {return locs_.size();}
  const double & operator[](const size_t ii) const {return locs_.at(ii);}

  const std::vector<double> & latitudes() const {return dummy_;}
  const std::vector<double> & longitudes() const {return locs_;}
  const std::vector<util::DateTime> & times() const {return times_;}

 private:
  void print(std::ostream & os) const;
  std::vector<double> locs_;
  std::vector<double> dummy_;  // for interfaces expecting a latitude coordinate
  std::vector<util::DateTime> times_;
};

}  // namespace lorenz95

#endif  // LORENZ95_LOCSL95_H_
