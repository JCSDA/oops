/*
 * (C) Copyright 2009-2016 ECMWF.
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
  typedef typename ObsAuxControl_::Parameters_ Parameters_;

  static const std::string classname() {return "oops::ObsAuxControl";}

  /// Constructor for specified ObsSpace \p os and \p params
  ObsAuxControl(const ObsSpace<OBS> & os, const Parameters_ & params);
  /// Creates ObsAuxControl with the same structure as \p other.
  /// Copies \p other if \p copy is true, otherwise creates zero ObsAuxControl
  explicit ObsAuxControl(const ObsAuxControl &, const bool copy = true);
  /// Destructor (defined explicitly for timing and tracing)
  ~ObsAuxControl();

  /// const Accessor
  const ObsAuxControl_ & obsauxcontrol() const {return *aux_;}
  /// Accessor
  ObsAuxControl_ & obsauxcontrol() {return *aux_;}

  /// Read this ObsAuxControl from file
  void read(const Parameters_ &);
  /// Write this ObsAuxControl out to file
  void write(const Parameters_ &) const;
  /// Norm (used in tests)
  double norm() const;

  /// Return required inputs variables from Model
  const Variables & requiredVars() const;
  /// Return required observations diagnostics
  const Variables & requiredHdiagnostics() const;

  /// Assign operator from other ObsAuxControl \p rhs
  ObsAuxControl & operator=(const ObsAuxControl & rhs);

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ObsAuxControl_> aux_;
};

// =============================================================================

template<typename OBS>
ObsAuxControl<OBS>::ObsAuxControl(const ObsSpace<OBS> & os,
                                    const Parameters_ & params) : aux_()
{
  Log::trace() << "ObsAuxControl<OBS>::ObsAuxControl starting" << std::endl;
  util::Timer timer(classname(), "ObsAuxControl");
  aux_.reset(new ObsAuxControl_(os.obsspace(), params));
  Log::trace() << "ObsAuxControl<OBS>::ObsAuxControl done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
ObsAuxControl<OBS>::ObsAuxControl(const ObsAuxControl & other, const bool copy) : aux_()
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
void ObsAuxControl<OBS>::read(const Parameters_ & params) {
  Log::trace() << "ObsAuxControl<OBS>::read starting" << std::endl;
  util::Timer timer(classname(), "read");
  aux_->read(params);
  Log::trace() << "ObsAuxControl<OBS>::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsAuxControl<OBS>::write(const Parameters_ & params) const {
  Log::trace() << "ObsAuxControl<OBS>::write starting" << std::endl;
  util::Timer timer(classname(), "write");
  aux_->write(params);
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
const Variables & ObsAuxControl<OBS>::requiredHdiagnostics() const {
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
