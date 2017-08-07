/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_GENERIC_LOCALIZATIONNICAS_H_
#define OOPS_GENERIC_LOCALIZATIONNICAS_H_

#include <sstream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "util/Logger.h"
#include "oops/interface/LocalizationBase.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// NICAS localization matrix.

template<typename MODEL>
class LocalizationNICAS : public LocalizationBase<MODEL> {
  typedef typename MODEL::Geometry              Geometry_;
  typedef typename MODEL::Increment             Increment_;

 public:
  LocalizationNICAS(const Geometry_ &, const eckit::Configuration &);
  ~LocalizationNICAS();

  void multiply(Increment_ &) const;

 private:
  void print(std::ostream &) const;
};

// =============================================================================

template<typename MODEL>
LocalizationNICAS<MODEL>::LocalizationNICAS(const Geometry_ &, const eckit::Configuration &)
{
  Log::trace() << "LocalizationNICAS:LocalizationNICAS constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
LocalizationNICAS<MODEL>::~LocalizationNICAS() {
  Log::trace() << "LocalizationNICAS:~LocalizationNICAS destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationNICAS<MODEL>::multiply(Increment_ & dx) const {
  Log::trace() << "LocalizationNICAS:multiply starting" << std::endl;

  UnstructuredGrid incr;
  dx.convert_to(incr);

  Log::info() << "LocalizationNICAS:multiply doing nothing" << std::endl;

  dx.convert_from(incr);

  Log::trace() << "LocalizationNICAS:multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationNICAS<MODEL>::print(std::ostream & os) const {
  os << "LocalizationNICAS<MODEL>::print not implemeted yet";
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_LOCALIZATIONNICAS_H_
