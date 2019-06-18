/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_BASE_VARIABLES_H_
#define OOPS_BASE_VARIABLES_H_

#include <ostream>
#include <string>
#include <vector>


#include "eckit/config/LocalConfiguration.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------

class Variables : public util::Printable {
 public:
  static const std::string classname() {return "oops::Variables";}

  Variables();
  explicit Variables(const eckit::Configuration &);
  explicit Variables(const std::vector<std::string> &, const std::string & conv = "");

  ~Variables();

  Variables(const Variables &);
  Variables & operator+=(const Variables &);

  size_t size() const {return vars_.size();}
  const std::string & operator[](const size_t kk) const {return vars_.at(kk);}

  bool has(const std::string &) const;
  size_t find(const std::string &) const;

  const std::vector<std::string> & variables() const {return vars_;}
  const std::vector<int> & channels() const {return channels_;}
  const eckit::Configuration & toFortran() const {return fconf_;}  // to be removed
  const eckit::Configuration * toFortranBetter() const {return &conf_;}

 private:
  void print(std::ostream &) const;
  void setConf();

  std::string convention_;
  std::vector<std::string> vars_;
  std::vector<int> channels_;        // channel indices
  eckit::LocalConfiguration conf_;
  eckit::LocalConfiguration fconf_;  // Until we can read vector of strings from fortran
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_VARIABLES_H_
