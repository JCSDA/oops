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
///        Another option is either construct the variable object with meta data
///        or add metadata to the variable object. The metadata keys typically are
///        variable names but don't need to be consistent with the vector of variable
///        names.
///
///        Most operators / methods do not affect/interact with the meta data.
///        The exceptions are:
///         += Variables;  where it updates the metadata and appends extra metadata
///                        from the right hand side Variables object.
///         == Variables;  where it compares the metadata for metadata keys that
///                        are internally consistent with the variable names.
///         .addMetaData(  that will either update a value or add a value within
///                        some metadata for a metadata (variable name) key.
class Variables : public util::Printable {
 public:
  static const std::string classname() {return "oops::Variables";}

  Variables();
  Variables(const eckit::Configuration &, const std::string &);
  explicit Variables(const std::vector<std::string> &, const std::string & conv = "");
  Variables(const std::vector<std::string> & vars, const std::vector<int> & channels);
  Variables(const eckit::Configuration &, const std::vector<std::string> & vars);

  Variables(const Variables &);
  Variables & operator+=(const Variables &);
  Variables & operator-=(const Variables &);
  Variables & operator-=(const std::string &);

  size_t size() const {return vars_.size();}
  const std::string & operator[](const size_t kk) const {return vars_.at(kk);}
  bool operator==(const Variables &) const;
  bool operator!=(const Variables &) const;
  bool operator<=(const Variables &) const;

  void addMetaData(const std::string & varname,
                   const std::string & keyname,
                   const int & keyvalue);

  bool has(const std::string &) const;

  size_t find(const std::string &) const;

  /// make this Variables an intersection between this Variables and other variables
  void intersection(const Variables & other);

  const std::vector<std::string> & variables() const {return vars_;}
  const std::vector<int> & channels() const {return channels_;}
  const eckit::Configuration & variablesMetaData() const {return varMetaData_;}

  void push_back(const std::string &);
  void sort();

  int getLevels(const std::string &) const;

 private:
  void print(std::ostream &) const;
  void setConf();
  /// returns sorted variable names
  std::vector<std::string> asCanonical() const;

  void getVariableSubKeyValue(const std::string & varname,
                              const std::string & keyname,
                              const eckit::Configuration & conf,
                              int & intvalue) const;

  void setVariableSubKeyValue(const std::string & varname,
                              const std::string & keyname,
                              const int & keyvalue,
                              eckit::LocalConfiguration & lconf);

  std::string convention_;
  std::vector<std::string> vars_;
  std::vector<int> channels_;        // channel indices
  eckit::LocalConfiguration varMetaData_;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_VARIABLES_H_
