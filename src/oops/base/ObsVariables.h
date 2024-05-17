/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "oops/base/VariablesBase.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Class to set an container of variable names and manipulate it.
///        One option on construction is to provide a vector of integers to specify
///        a satellite channel list selection. These integers are convolved with
///        with original vector of variable names to create a new vector of variable
///        names (original variable name + "_" + channel number).
class ObsVariables : public VariablesBase {
 public:
  static const std::string classname() {return "oops::ObsVariables";}

  ObsVariables() = default;
  ObsVariables(const eckit::Configuration &, const std::string &);
  explicit ObsVariables(const std::vector<std::string> &);
  ObsVariables(const std::vector<std::string> & vars, const std::vector<int> & channels);

  ObsVariables & operator+=(const ObsVariables &);
  ObsVariables & operator-=(const std::string &);

  bool operator==(const ObsVariables &) const;
  bool operator!=(const ObsVariables &) const;

  /// make this ObsVariables an intersection between this ObsVariables and other ObsVariables
  void intersection(const ObsVariables & other);

  const std::vector<int> & channels() const {return channels_;}

  void sort();

 private:
  std::vector<int> channels_;        // channel indices
};

// -----------------------------------------------------------------------------

}  // namespace oops
