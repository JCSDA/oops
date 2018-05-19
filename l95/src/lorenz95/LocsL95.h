/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_LOCSL95_H_
#define LORENZ95_LOCSL95_H_

#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace lorenz95 {

/// LocsL95 class to handle locations for L95 model.
class LocsL95 : public util::Printable,
                private util::ObjectCounter<LocsL95> {
 public:
  static const std::string classname() {return "lorenz95::LocsL95";}

  LocsL95(const std::vector<int> &, const std::vector<double> &);
  explicit LocsL95(const eckit::Configuration &);
  ~LocsL95() {}

  size_t size() const {return locs_.size();}
  const double & operator[](const size_t ii) const {return locs_.at(ii);}
  const int & globalIndex(const size_t ii) const {return indx_.at(ii);}

 private:
  void print(std::ostream & os) const;
  std::vector<int> indx_;
  std::vector<double> locs_;
};

}  // namespace lorenz95

#endif  // LORENZ95_LOCSL95_H_
