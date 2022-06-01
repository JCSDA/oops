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

#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

class Variables : public util::Printable {
 public:
  static const std::string classname() {return "oops::Variables";}

  Variables();
  Variables(const eckit::Configuration &, const std::string &);
  explicit Variables(const std::vector<std::string> &, const std::string & conv = "");
  Variables(const std::vector<std::string> & vars, const std::vector<int> & channels);

  Variables(const Variables &);
  Variables & operator+=(const Variables &);
  Variables & operator-=(const Variables &);
  Variables & operator-=(const std::string &);

  size_t size() const {return vars_.size();}
  const std::string & operator[](const size_t kk) const {return vars_.at(kk);}
  bool operator==(const Variables &) const;
  bool operator!=(const Variables &) const;
  bool operator<=(const Variables &) const;

  bool has(const std::string &) const;
  size_t find(const std::string &) const;

  /// make this Variables an intersection between this Variables and \p other variables
  void intersection(const Variables & other);

  const std::vector<std::string> & variables() const {return vars_;}
  const std::vector<int> & channels() const {return channels_;}
  void push_back(const std::string &);
  void sort();

 private:
  void print(std::ostream &) const;
  void setConf();
  /// returns sorted variable names
  std::vector<std::string> asCanonical() const;

  std::string convention_;
  std::vector<std::string> vars_;
  std::vector<int> channels_;        // channel indices
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_VARIABLES_H_
