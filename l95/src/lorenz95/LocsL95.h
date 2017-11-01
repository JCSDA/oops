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

#include "lorenz95/ObsTable.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace lorenz95 {
  class ObservationL95;

/// LocsL95 class to handle locations for L95 model.

class LocsL95 : public util::Printable,
                private util::ObjectCounter<LocsL95> {
 public:
  static const std::string classname() {return "lorenz95::LocsL95";}

  LocsL95(const ObsTable & ot,
          const util::DateTime & t1, const util::DateTime & t2) {
    locs_ = ot.locations(t1, t2);
  }
  ~LocsL95() {}

  int nobs() const {return locs_.size();}

  const double & operator[](const int ii) const {return locs_[ii];}

 private:
  void print(std::ostream & os) const {
    os << "LocsL95::print not implemented";
  }
  std::vector<double> locs_;
};

}  // namespace lorenz95

#endif  // LORENZ95_LOCSL95_H_
