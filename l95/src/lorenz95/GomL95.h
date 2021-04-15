/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_GOML95_H_
#define LORENZ95_GOML95_H_

#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace oops {
  class Variables;
}

namespace lorenz95 {
  class LocsL95;
  class ObsTable;

/// GomL95 class to handle locations for L95 model.

class GomL95 : public util::Printable,
               private util::ObjectCounter<GomL95> {
 public:
  static const std::string classname() {return "lorenz95::GomL95";}

  GomL95(const LocsL95 &, const oops::Variables &);
  GomL95(const eckit::Configuration &, const ObsTable &,
         const oops::Variables &);

  void zero();
  void random();
  double rms() const;
  double normalizedrms(const GomL95 &) const;
  GomL95 & operator*=(const double &);
  GomL95 & operator+=(const GomL95 &);
  GomL95 & operator-=(const GomL95 &);
  GomL95 & operator*=(const GomL95 &);
  double dot_product_with(const GomL95 &) const;
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  void print(std::ostream &) const;

  size_t size() const {return size_;}
  const double & operator[](const int ii) const {return locval_[ii];}
  double & operator[](const int ii) {return locval_[ii];}

 private:
  size_t size_;
  std::vector<double> locval_;
};

}  // namespace lorenz95

#endif  // LORENZ95_GOML95_H_
