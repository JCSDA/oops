/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_MODELAUXCOVARIANCE_H_
#define OOPS_INTERFACE_MODELAUXCOVARIANCE_H_

#include <iostream>
#include <memory>
#include <string>
#include <type_traits>

#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/Geometry.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Auxiliary Error Covariance related to model, not used at the moment.
/// \details
/// This class calls the model's implementation of ModelAuxCovariance.
// -----------------------------------------------------------------------------

template <typename MODEL>
class ModelAuxCovariance : public util::Printable,
                           private boost::noncopyable,
                           private util::ObjectCounter<ModelAuxCovariance<MODEL> > {
  typedef typename MODEL::ModelAuxCovariance    ModelAuxCovariance_;
  typedef Geometry<MODEL>             Geometry_;
  typedef ModelAuxControl<MODEL>      ModelAuxControl_;
  typedef ModelAuxIncrement<MODEL>    ModelAuxIncrement_;

 public:
  static const std::string classname() {return "oops::ModelAuxCovariance";}

  ModelAuxCovariance(const eckit::Configuration &, const Geometry_ &);
  /// Destructor (defined explicitly for timing and tracing)
  ~ModelAuxCovariance();

  /// linearize operator
  void linearize(const ModelAuxControl_ &, const Geometry_ &);
  /// Sets the second parameter to the first multiplied by the covariance matrix.
  void multiply(const ModelAuxIncrement_ &, ModelAuxIncrement_ &) const;
  /// Sets the second parameter to the first multiplied by the inverse covariance matrix.
  void inverseMultiply(const ModelAuxIncrement_ &, ModelAuxIncrement_ &) const;
  /// randomize the values in the ModelAuxIncrement
  void randomize(ModelAuxIncrement_ &) const;

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ModelAuxCovariance_> cov_;
};

// =============================================================================

template<typename MODEL>
ModelAuxCovariance<MODEL>::ModelAuxCovariance(const eckit::Configuration & conf,
                                              const Geometry_ & resol) : cov_()
{
  Log::trace() << "ModelAuxCovariance<MODEL>::ModelAuxCovariance starting" << std::endl;
  util::Timer timer(classname(), "ModelAuxCovariance");
  cov_.reset(new ModelAuxCovariance_(conf, resol.geometry()));
  Log::trace() << "ModelAuxCovariance<MODEL>::ModelAuxCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ModelAuxCovariance<MODEL>::~ModelAuxCovariance() {
  Log::trace() << "ModelAuxCovariance<MODEL>::~ModelAuxCovariance starting" << std::endl;
  util::Timer timer(classname(), "~ModelAuxCovariance");
  cov_.reset();
  Log::trace() << "ModelAuxCovariance<MODEL>::~ModelAuxCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ModelAuxCovariance<MODEL>::linearize(const ModelAuxControl_ & xx, const Geometry_ & resol) {
  Log::trace() << "ModelAuxCovariance<MODEL>::linearize starting" << std::endl;
  util::Timer timer(classname(), "linearize");
  cov_->linearize(xx.modelauxcontrol(), resol.geometry());
  Log::trace() << "ModelAuxCovariance<MODEL>::linearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ModelAuxCovariance<MODEL>::multiply(const ModelAuxIncrement_ & dx1,
                                         ModelAuxIncrement_ & dx2) const {
  Log::trace() << "ModelAuxCovariance<MODEL>::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  cov_->multiply(dx1.modelauxincrement(), dx2.modelauxincrement());
  Log::trace() << "ModelAuxCovariance<MODEL>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ModelAuxCovariance<MODEL>::inverseMultiply(const ModelAuxIncrement_ & dx1,
                                                ModelAuxIncrement_ & dx2) const {
  Log::trace() << "ModelAuxCovariance<MODEL>::inverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "inverseMultiply");
  cov_->inverseMultiply(dx1.modelauxincrement(), dx2.modelauxincrement());
  Log::trace() << "ModelAuxCovariance<MODEL>::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ModelAuxCovariance<MODEL>::randomize(ModelAuxIncrement_ & dx) const {
  Log::trace() << "ModelAuxCovariance<MODEL>::randomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");
  cov_->randomize(dx.modelauxincrement());
  Log::trace() << "ModelAuxCovariance<MODEL>::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ModelAuxCovariance<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ModelAuxCovariance<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *cov_;
  Log::trace() << "ModelAuxCovariance<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_MODELAUXCOVARIANCE_H_
