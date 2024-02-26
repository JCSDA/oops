/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_WRITEPARAMETERSBASE_H_
#define OOPS_BASE_WRITEPARAMETERSBASE_H_

#include <string>

#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Base class of classes storing options controlling the process of writing a state or
/// increment to a file.
class WriteParametersBase : public Parameters {
  OOPS_ABSTRACT_PARAMETERS(WriteParametersBase, Parameters)

 public:
  /// Set the ensemble member index to \p mem.
  void setMember(const int &);

  /// Set the validity date/time to \p d.
  void setDate(const util::DateTime &);

  /// Set the Krylov solver iteration to \p it.
  void setIteration(const int &);

  /// Type
  OptionalParameter<std::string> type{"type", this};

  /// Experiment id
  OptionalParameter<std::string> exp{"exp", this};

  /// Ensemble member index.
  OptionalParameter<int> member{"member", this};

  /// Ensemble member pattern.
  OptionalParameter<std::string> memberPattern{"member pattern", this};

  /// Validity time.
  OptionalParameter<util::DateTime> date{"date", this};

  /// Krylov solver iteration.
  OptionalParameter<int> iteration{"iteration", this};

  /// Prefix for filenames
  OptionalParameter<std::string> prefix{"prefix", this};

  /// IO format for the date in the filenames
  Parameter<bool> dateCols{"date colons", true, this};
};

// -----------------------------------------------------------------------------

/// \brief A subclass of WriteParametersBase storing the values of all options in an
/// eckit::LocalConfiguration object.
///
/// This object can be accessed by calling the `value()` method of the `config` member variable.
///
/// Values loaded into a ConfigurationParameter aren't validated in any way; parts of the code
/// using GenericWriteParameters should therefore ideally be refactored, replacing this class with
/// a dedicated subclass of Parameters storing each parameter in a separate
/// (Optional/Required)Parameter object.
class GenericWriteParameters : public WriteParametersBase {
  OOPS_CONCRETE_PARAMETERS(GenericWriteParameters, WriteParametersBase)

 public:
  ConfigurationParameter config{this};
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_WRITEPARAMETERSBASE_H_
