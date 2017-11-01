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

#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace util {
  class DateTime;
}

namespace lorenz95 {
  class ObsTable;
  class NoVariables;
  class Resolution;

/// GomL95 class to handle locations for L95 model.

class GomL95 : public util::Printable,
               private util::ObjectCounter<GomL95> {
 public:
  static const std::string classname() {return "lorenz95::GomL95";}

  GomL95(const ObsTable &, const NoVariables &, const util::DateTime &, const util::DateTime &, const Resolution &);
  ~GomL95();

  double dot_product_with(const GomL95 &) const;
  void zero();

  const double & operator[](const int ii) const {return locval_[ii];}
  double & operator[](const int ii) {return locval_[ii];}
  int getindx(const int il) const {return iobs_[il];}
  int nobs() const {return size_;}
  int & current() const {return current_;}

 private:
  void print(std::ostream &) const;
  int size_;
  std::vector<int> iobs_;
  std::vector<double> locval_;
  mutable int current_;
};

}  // namespace lorenz95

#endif  // LORENZ95_GOML95_H_
