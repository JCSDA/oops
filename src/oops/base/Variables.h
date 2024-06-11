/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variable.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Class to set an container of Variable objects.
///        If constructed only from a vector of strings, the metadata for each Variable is
///        defaulted, and no vertical levels are set.
///        If constructed from a vector of vector of strings and a Configuration object, the
///        levels can be passed in the Configuration object.
///        This class has been designed to not allow duplicate variables, so the push_back method
///        will not add a variable if it already exists in the container.
class Variables : public util::Printable {
 public:
  static const std::string classname() {return "oops::Variables";}

  Variables() = default;
  Variables(const eckit::Configuration &, const std::string &);
  explicit Variables(const std::vector<std::string> &);
  explicit Variables(const std::vector<Variable> &);
  Variables(const eckit::Configuration &, const std::vector<std::string> & vars);

  size_t size() const {return vars_.size();}
  const oops::Variable & operator[](const size_t kk) const {return vars_.at(kk);}
  oops::Variable & operator[](const size_t kk) {return vars_.at(kk);}
  const oops::Variable & operator[](const std::string &) const;
  oops::Variable & operator[](const std::string &);
  // TODO(AS): this method needs to be removed.
  const std::vector<std::string> variables() const;

  bool has(const Variable &) const;
  size_t find(const Variable &) const;
  bool has(const std::string &) const;
  size_t find(const std::string &) const;
  void push_back(const Variable &);
  void push_back(const std::string &);

  Variables & operator+=(const Variables &);
  Variables & operator-=(const Variables &);
  Variables & operator-=(const Variable &);

  bool operator==(const Variables &) const;
  bool operator!=(const Variables &) const;
  bool operator<=(const Variables &) const;

  auto begin()  const { return vars_.begin(); }
  auto begin()        { return vars_.begin(); }
  auto cbegin() const { return vars_.cbegin(); }
  auto end()  const { return vars_.end(); }
  auto end()        { return vars_.end(); }
  auto cend() const { return vars_.cend(); }

  /// make this Variables an intersection between this Variables and other variables
  void intersection(const Variables & other);

  void sort();

 private:
  void print(std::ostream &) const override;
  /// returns sorted by variable names
  std::vector<Variable> asCanonical() const;

  /// Data
  std::vector<Variable> vars_;
};

// -----------------------------------------------------------------------------

}  // namespace oops
