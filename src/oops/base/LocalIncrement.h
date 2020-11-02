/*
 * (C) Copyright 2018-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_LOCALINCREMENT_H_
#define OOPS_BASE_LOCALINCREMENT_H_

#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

namespace oops {

class LocalIncrement: public util::Printable {
 public:
  LocalIncrement(const oops::Variables vars,
                 std::vector<double> vals,
                 std::vector<int> varlens)
       : vars_(vars), vals_(vals), varlens_(varlens) {}

  const oops::Variables & getVars() const {return vars_;}
  const std::vector<double> & getVals() const {return vals_;}
  void setVals(std::vector<double> &);

  /// Linear algebra operators
  LocalIncrement & operator*=(const std::vector<double> &);

 private:
  void print(std::ostream & os) const {
    os << "LocalIncrement, size: " << vals_.size() << ", first element: "
       << vals_[0] << std::endl; }
  const oops::Variables vars_;      // variables in the object
  std::vector<double> vals_;        // data in flat array
  const std::vector<int> varlens_;  // vector containing nlevs for each variable
};

}  // namespace oops

#endif  // OOPS_BASE_LOCALINCREMENT_H_
