/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Base class for a container of variable names.
class VariablesBase : public util::Printable {
 public:
  VariablesBase() = default;
  VariablesBase(const eckit::Configuration &, const std::string &);
  explicit VariablesBase(const std::vector<std::string> &);

  size_t size() const {return vars_.size();}
  const std::string & operator[](const size_t kk) const {return vars_.at(kk);}

  bool has(const std::string &) const;
  size_t find(const std::string &) const;

  const std::vector<std::string> & variables() const {return vars_;}

  void push_back(const std::string &);

 private:
  void print(std::ostream &) const override;

 protected:
  /// returns sorted variable names
  std::vector<std::string> asCanonical() const;
  std::vector<std::string> vars_;
};

// -----------------------------------------------------------------------------

}  // namespace oops
