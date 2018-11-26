/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_OBSERRORCOVARIANCE_H_
#define OOPS_INTERFACE_OBSERRORCOVARIANCE_H_

#include <string>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/interface/ObsErrorBase.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// Observation error covariance matrix
/*!
 *  This class provides the operations associated with the observation
 *  error covariance matrix. It wraps the actual observation error covariance
 *  which can be a model specific one or a generic one.
 *  The interface for the observation error comprises two levels (ObsErrorCovariance
 *  and ObsErrorBase) because we want run time polymorphism.
 *  The ObsErrorCovariance does conversion of arguments to templated ObsVector and
 *  the tracing and timing. The ObsErrorBase does the conversion to model specific
 *  ObsVector.
 */

template <typename MODEL>
class ObsErrorCovariance : public util::Printable,
                           private util::ObjectCounter<ObsErrorCovariance<MODEL> >,
                           private boost::noncopyable {
  typedef ObsErrorBase<MODEL>        ObsErrorBase_;
  typedef ObservationSpace<MODEL>    ObsSpace_;
  typedef ObsVector<MODEL>           ObsVector_;

 public:
  static const std::string classname() {return "oops::ObsErrorCovariance";}

  ObsErrorCovariance(const eckit::Configuration &, const ObsSpace_ &, const Variables &);
  ~ObsErrorCovariance();

/// Update after QC or other obs filters
  void update();

/// Multiply a Departure by \f$R\f$ and \f$R^{-1}\f$
  void multiply(ObsVector_ &) const;
  void inverseMultiply(ObsVector_ &) const;

/// Generate random perturbation
  void randomize(ObsVector_ &) const;

/// Get mean error for Jo table
  double getRMSE() const;

 private:
  void print(std::ostream &) const;
  boost::scoped_ptr<ObsErrorBase_> covar_;
};

// ====================================================================================

template <typename MODEL>
ObsErrorCovariance<MODEL>::ObsErrorCovariance(const eckit::Configuration & conf,
                                              const ObsSpace_ & obsdb,
                                              const Variables & obsvar) : covar_() {
  Log::trace() << "ObsErrorCovariance<MODEL>::ObsErrorCovariance starting" << std::endl;
  util::Timer timer(classname(), "ObsErrorCovariance");
  covar_.reset(ObsErrorFactory<MODEL>::create(conf, obsdb, obsvar));
  Log::trace() << "ObsErrorCovariance<MODEL>::ObsErrorCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsErrorCovariance<MODEL>::~ObsErrorCovariance() {
  Log::trace() << "ObsErrorCovariance<MODEL>::~ObsErrorCovariance starting" << std::endl;
  util::Timer timer(classname(), "~ObsErrorCovariance");
  covar_.reset();
  Log::trace() << "ObsErrorCovariance<MODEL>::~ObsErrorCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObsErrorCovariance<MODEL>::update() {
  Log::trace() << "ObsErrorCovariance<MODEL>::update starting" << std::endl;
  util::Timer timer(classname(), "update");
  covar_->update();
  Log::trace() << "ObsErrorCovariance<MODEL>::update done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObsErrorCovariance<MODEL>::multiply(ObsVector_ & dy) const {
  Log::trace() << "ObsErrorCovariance<MODEL>::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  covar_->multiply(dy.obsvector());
  Log::trace() << "ObsErrorCovariance<MODEL>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObsErrorCovariance<MODEL>::inverseMultiply(ObsVector_ & dy) const {
  Log::trace() << "ObsErrorCovariance<MODEL>::inverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "inverseMultiply");
  covar_->inverseMultiply(dy.obsvector());
  Log::trace() << "ObsErrorCovariance<MODEL>::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObsErrorCovariance<MODEL>::randomize(ObsVector_ & dy) const {
  Log::trace() << "ObsErrorCovariance<MODEL>::randomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");
  covar_->randomize(dy.obsvector());
  Log::trace() << "ObsErrorCovariance<MODEL>::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
double ObsErrorCovariance<MODEL>::getRMSE() const {
  Log::trace() << "ObsErrorCovariance<MODEL>::getRMSE starting" << std::endl;
  util::Timer timer(classname(), "getRMSE");
  double zz = covar_->getRMSE();
  Log::trace() << "ObsErrorCovariance<MODEL>::getRMSE done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsErrorCovariance<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ObsErrorCovariance<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *covar_;
  Log::trace() << "ObsErrorCovariance<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSERRORCOVARIANCE_H_
