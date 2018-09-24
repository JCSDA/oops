/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_GRIDPOINT_H_
#define OOPS_BASE_GRIDPOINT_H_

#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

namespace oops {

class GridPoint: public util::Printable {
 public:
  GridPoint(const oops::Variables vars, std::vector<double> vals, std::vector<int> varlens)
       : vars_(vars), vals_(vals), varlens_(varlens) {}
  ~GridPoint() {}

  const std::vector<double> getVals() {return vals_;}
 private:
  void print(std::ostream & os) const {
    os << "GridPoint, size: " << vals_.size() << ", first element: " << vals_[0] << std::endl; }
  const oops::Variables vars_;
  const std::vector<double> vals_;     // data in flat array
  const std::vector<int>    varlens_;  // vector containing nlevs for each variable
};

}  // namespace oops

#endif  // OOPS_BASE_GRIDPOINT_H_
