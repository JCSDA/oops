/*
 * (C) Copyright 2018  UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_OBSDIAGNOSTICS_H_
#define OOPS_INTERFACE_OBSDIAGNOSTICS_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/base/ObsVariables.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

template <typename OBS> class Locations;

// -----------------------------------------------------------------------------

template <typename OBS>
class ObsDiagnostics : public util::Printable,
                       private boost::noncopyable,
                       private util::ObjectCounter<ObsDiagnostics<OBS> > {
  typedef typename OBS::ObsDiagnostics   ObsDiags_;
  typedef ObsSpace<OBS>                  ObsSpace_;
  typedef Locations<OBS>                 Locations_;

 public:
  static const std::string classname() {return "oops::ObsDiagnostics";}

  ObsDiagnostics(const ObsSpace_ &, const Locations_ &, const ObsVariables &);
  // ctor used in the test for ObsFilters (not implemented in toy models)
  ObsDiagnostics(const eckit::Configuration &, const ObsSpace_ &, const ObsVariables &);

  ~ObsDiagnostics();

/// Interfacing
  ObsDiags_ & obsdiagnostics() {return *diags_;}
  const ObsDiags_ & obsdiagnostics() const {return *diags_;}

// I/O
  void save(const std::string &) const;

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ObsDiags_> diags_;
};

// -----------------------------------------------------------------------------
template <typename OBS>
ObsDiagnostics<OBS>::ObsDiagnostics(const ObsSpace_ & os,
                                    const Locations_ & locations,
                                    const ObsVariables & vars) : diags_()
{
  Log::trace() << "ObsDiagnostics<OBS>::ObsDiagnostics starting" << std::endl;
  util::Timer timer(classname(), "ObsDiagnostics");
  diags_.reset(new ObsDiags_(os.obsspace(), locations, vars));
  Log::trace() << "ObsDiagnostics<OBS>::ObsDiagnostics done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS>
ObsDiagnostics<OBS>::ObsDiagnostics(const eckit::Configuration & config, const ObsSpace_ & os,
                                    const ObsVariables & vars) : diags_()
{
  Log::trace() << "ObsDiagnostics<OBS>::ObsDiagnostics starting" << std::endl;
  util::Timer timer(classname(), "ObsDiagnostics");
  diags_.reset(new ObsDiags_(config, os.obsspace(), vars));
  Log::trace() << "ObsDiagnostics<OBS>::ObsDiagnostics done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS>
ObsDiagnostics<OBS>::~ObsDiagnostics() {
  Log::trace() << "ObsDiagnostics<OBS>::~ObsDiagnostics starting" << std::endl;
  util::Timer timer(classname(), "~ObsDiagnostics");
  diags_.reset();
  Log::trace() << "ObsDiagnostics<OBS>::~ObsDiagnostics done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS>
void ObsDiagnostics<OBS>::save(const std::string & name) const {
  Log::trace() << "ObsDiagnostics<OBS, DATATYPE>::save starting " << name << std::endl;
  util::Timer timer(classname(), "save");
  diags_->save(name);
  Log::trace() << "ObsDiagnostics<OBS, DATATYPE>::save done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS>
void ObsDiagnostics<OBS>::print(std::ostream & os) const {
  Log::trace() << "ObsDiagnostics<OBS>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *diags_;
  Log::trace() << "ObsDiagnostics<OBS>::print done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSDIAGNOSTICS_H_
