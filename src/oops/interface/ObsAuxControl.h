/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Crown Copyright 2023, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_OBSAUXCONTROL_H_
#define OOPS_INTERFACE_OBSAUXCONTROL_H_

#include <iostream>
#include <memory>
#include <string>

#include "oops/interface/ObsSpace.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {
  class ObsVariables;
  class Variables;

// -----------------------------------------------------------------------------
/// \brief Auxiliary state related to observations, templated on <OBS>
/// \details
/// This is currently only used for bias correction coefficients, but can be used for other cases.
/// This class calls the <OBS> implementation of ObsAuxControl.
// -----------------------------------------------------------------------------

template <typename OBS>
class ObsAuxControl : public util::Printable,
                      private util::ObjectCounter<ObsAuxControl<OBS> > {
  typedef typename OBS::ObsAuxControl          ObsAuxControl_;

 public:
  static const std::string classname() {return "oops::ObsAuxControl";}

  /// Constructor for specified ObsSpace \p os and \p params
  ObsAuxControl(const ObsSpace<OBS> & os, const eckit::Configuration &);
  /// Creates ObsAuxControl with the same structure as \p other.
  /// Copies \p other if \p copy is true, otherwise creates zero ObsAuxControl
  explicit ObsAuxControl(const ObsAuxControl &, const bool copy = true);
  /// Destructor (defined explicitly for timing and tracing)
  ~ObsAuxControl();

  /// const Accessor
  const ObsAuxControl_ & obsauxcontrol() const {return *aux_;}
  /// Accessor
  ObsAuxControl_ & obsauxcontrol() {return *aux_;}

  const eckit::Configuration & config() const {return config_;}
  const ObsSpace<OBS> & obspace() const {return obspace_;}

  /// Read this ObsAuxControl from file
  void read(const eckit::Configuration &);
  /// Write this ObsAuxControl out to file
  void write(const eckit::Configuration &) const;
  /// Norm (used in tests)
  double norm() const;

  /// Return required inputs variables from Model
  const Variables & requiredVars() const;
  /// Return required observations diagnostics
  const ObsVariables & requiredHdiagnostics() const;

  /// Assign operator from other ObsAuxControl \p rhs
  ObsAuxControl & operator=(const ObsAuxControl & rhs);

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ObsAuxControl_> aux_;
  const eckit::LocalConfiguration config_;
  const ObsSpace<OBS> & obspace_;
};

// =============================================================================

template<typename OBS>
ObsAuxControl<OBS>::ObsAuxControl(const ObsSpace<OBS> & os, const eckit::Configuration & config)
  : aux_(), config_(config), obspace_(os)
{
  Log::trace() << "ObsAuxControl<OBS>::ObsAuxControl starting" << std::endl;
  util::Timer timer(classname(), "ObsAuxControl");
  aux_.reset(new ObsAuxControl_(os.obsspace(), config));
  Log::trace() << "ObsAuxControl<OBS>::ObsAuxControl done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
ObsAuxControl<OBS>::ObsAuxControl(const ObsAuxControl & other, const bool copy)
  : aux_(), config_(other.config_), obspace_(other.obspace_)
{
  Log::trace() << "ObsAuxControl<OBS>::ObsAuxControl copy starting" << std::endl;
  util::Timer timer(classname(), "ObsAuxControl");
  aux_.reset(new ObsAuxControl_(*other.aux_, copy));
  Log::trace() << "ObsAuxControl<OBS>::ObsAuxControl copy done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
ObsAuxControl<OBS>::~ObsAuxControl() {
  Log::trace() << "ObsAuxControl<OBS>::~ObsAuxControl starting" << std::endl;
  util::Timer timer(classname(), "~ObsAuxControl");
  aux_.reset();
  Log::trace() << "ObsAuxControl<OBS>::~ObsAuxControl done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsAuxControl<OBS>::read(const eckit::Configuration & config) {
  Log::trace() << "ObsAuxControl<OBS>::read starting" << std::endl;
  util::Timer timer(classname(), "read");
  aux_->read(config);
  Log::trace() << "ObsAuxControl<OBS>::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsAuxControl<OBS>::write(const eckit::Configuration & config) const {
  Log::trace() << "ObsAuxControl<OBS>::write starting" << std::endl;
  util::Timer timer(classname(), "write");
  aux_->write(config);
  Log::trace() << "ObsAuxControl<OBS>::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
double ObsAuxControl<OBS>::norm() const {
  Log::trace() << "ObsAuxControl<OBS>::norm starting" << std::endl;
  util::Timer timer(classname(), "norm");
  double zz = aux_->norm();
  Log::trace() << "ObsAuxControl<OBS>::norm done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename OBS>
const Variables & ObsAuxControl<OBS>::requiredVars() const {
  Log::trace() << "ObsAuxControl<OBS>::requiredVars starting" << std::endl;
  util::Timer timer(classname(), "requiredVars");
  Log::trace() << "ObsAuxControl<OBS>::requiredVars done" << std::endl;
  return aux_->requiredVars();
}

// -----------------------------------------------------------------------------

template<typename OBS>
const ObsVariables & ObsAuxControl<OBS>::requiredHdiagnostics() const {
  Log::trace() << "ObsAuxControl<OBS>::requiredHdiagnostics starting" << std::endl;
  util::Timer timer(classname(), "requiredHdiagnostics");
  Log::trace() << "ObsAuxControl<OBS>::requiredHdiagnostics done" << std::endl;
  return aux_->requiredHdiagnostics();
}

// -----------------------------------------------------------------------------
template<typename OBS>
ObsAuxControl<OBS> & ObsAuxControl<OBS>::operator=(const ObsAuxControl & rhs) {
  Log::trace() << "ObsAuxControl<OBS>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  *aux_ = *rhs.aux_;
  Log::trace() << "ObsAuxControl<OBS>::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsAuxControl<OBS>::print(std::ostream & os) const {
  Log::trace() << "ObsAuxControl<OBS>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *aux_;
  Log::trace() << "ObsAuxControl<OBS>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSAUXCONTROL_H_
