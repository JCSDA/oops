/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_BASE_VARIABLES_H_
#define OOPS_BASE_VARIABLES_H_

#include <ostream>
#include <string>
#include <vector>

#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

class Variables : public util::Printable {
 public:
  static const std::string classname() {return "oops::Variables";}

  explicit Variables(const eckit::Configuration &);
  Variables(const std::vector<std::string> &, const std::string & conv = "");

  ~Variables();

  Variables(const Variables &);

  const int * toFortran() const {return &fvars_[0];}

 private:
  void print(std::ostream &) const;
  void setF90();

  std::string convention_;
  std::vector<std::string> vars_;
  std::vector<int> fvars_;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_VARIABLES_H_
