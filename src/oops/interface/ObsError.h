/*
 * (C) Crown copyright 2021, Met Office
 * (C) Copyright 2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_OBSERROR_H_
#define OOPS_INTERFACE_OBSERROR_H_

#include <memory>
#include <string>

#include "oops/interface/ObsSpace.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Wrapper class for OBS-specific implementations of the ObsError interface.
template <typename OBS>
class ObsError : public util::Printable,
                 private util::ObjectCounter<ObsError<OBS> > {
  typedef typename OBS::ObsError  ObsError_;

 public:
  // the obsSpace passed to this constructor is an oops wrapper around an obs-specific
  // obs space. Everything in oops/interface is a wrapper around a obs/model specific implementation
  ObsError(const eckit::Configuration & conf, const ObsSpace<OBS> & os) : timeComm_(os.timeComm()) {
    Log::trace() << "ObsError<OBS>::ObsError starting" << std::endl;
    // making an obs-specific ObsError class, need to get the
    // obs-specific obs-space from the oops obsspace wrapper (os).
    error_ = std::make_unique<ObsError_> (conf, os.obsspace());

    Log::trace() << "ObsError<OBS>::ObsError done" << std::endl;
  }

  ~ObsError();
  ObsError(const ObsError &) = delete;
  ObsError(ObsError &&) = default;
  ObsError& operator=(const ObsError &) = delete;
  ObsError& operator=(ObsError &&) = default;

  static const std::string classname() {return "oops::ObsError";}

  // Wrapper methods to access the obs-specific implementations in the wrapped error_ object

  // dy is an oops interface wrapper around an obs-specific obs vector
  // implementation. dy.obsvector() returns the obs-specific implementation
  // of an obsVector (eg. ioda::ObsVector or lorenz95::ObsVec1D)
  void multiply(oops::ObsVector<OBS> &dy) const {
    Log::trace() << "ObsError<OBS>::multiply starting" << std::endl;
    util::Timer timer(classname(), "multiply");
    error_->multiply(dy.obsvector());
    Log::trace() << "ObsError<OBS>::multiply done" << std::endl;
  }

  void inverseMultiply(oops::ObsVector<OBS> &dy) const {
    Log::trace() << "ObsError<OBS>::inverseMultiply starting" << std::endl;
    util::Timer timer(classname(), "inverseMultiply");
    error_->inverseMultiply(dy.obsvector());
    Log::trace() << "ObsError<OBS>::inverseMultiply done" << std::endl;
  }

  void randomize(oops::ObsVector<OBS> &dy) const {
    Log::trace() << "ObsError<OBS>::randomize starting" << std::endl;
    util::Timer timer(classname(), "randomize");
    error_->randomize(dy.obsvector());
    Log::trace() << "ObsError<OBS>::randomize done" << std::endl;
  }

  void save(const std::string & name) const {
    Log::trace() << "ObsError<OBS>::save starting" << std::endl;
    util::Timer timer(classname(), "save");
    error_->save(name);
    Log::trace() << "ObsError<OBS>::save done" << std::endl;
  }

  double getRMSE() const {
    Log::trace() << "ObsError<OBS>::getRMSE starting" << std::endl;
    util::Timer timer(classname(), "getRMSE");
    // get rmse from wrapped oops obsvector which knows how
    // to do time communication across 4denvar sub windows
    return this->obserrors().rms();
  }

  oops::ObsVector<OBS> obserrors() const {
    Log::trace() << "ObsError<OBS>::obserrors starting" << std::endl;
    util::Timer timer(classname(), "obserrors");
    return oops::ObsVector<OBS>(error_->getObsErrors(), timeComm_);
  }

  void update(const oops::ObsVector<OBS> &dy) {
    Log::trace() << "ObsError<OBS>::update starting" << std::endl;
    util::Timer timer(classname(), "update");
    error_->update(dy.obsvector());
    Log::trace() << "ObsError<OBS>::update done" << std::endl;
  }

  oops::ObsVector<OBS> inverseVariance() const {
    Log::trace() << "ObsError<OBS>::inverseVariance starting" << std::endl;
    util::Timer timer(classname(), "inverseVariance");
    return oops::ObsVector<OBS>(error_->getInverseVariance(), timeComm_);
  }

  void localize(oops::ObsVector<OBS> &locvector) const {
    Log::trace() << "ObsError<OBS>::localize starting" << std::endl;
    util::Timer timer(classname(), "localize");
    error_->localize(locvector.obsvector());
    Log::trace() << "ObsError<OBS>::localize done" << std::endl;
  }

  int localDim() const {
    Log::trace() << "ObsError<OBS>::localDim starting" << std::endl;
    util::Timer timer(classname(), "localDim");
    return error_->localDim();
  }

  Eigen::MatrixXf localInverseMultiply(const Eigen::MatrixXf &zz) const {
    Log::trace() << "ObsError<OBS>::localInverseMultiply starting" << std::endl;
    util::Timer timer(classname(), "localInverseMultiply");
    return error_->localInverseMultiply(zz);
  }

  Eigen::VectorXd local_invVarR() const {
    Log::trace() << "ObsError<OBS>::local_invVarR starting" << std::endl;
    util::Timer timer(classname(), "local_invVarR");
    return error_->local_invVarR();
  }

 private:
  void print(std::ostream &) const;

  std::unique_ptr<ObsError_> error_;
  const eckit::mpi::Comm & timeComm_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsError<OBS>::~ObsError() {
  Log::trace() << "ObsError<OBS>::~ObsError starting" << std::endl;
  util::Timer timer(classname(), "~ObsError");
  error_.reset();
  Log::trace() << "ObsError<OBS>::~ObsError done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsError<OBS>::print(std::ostream & os) const {
  Log::trace() << "ObsError<OBS>::print starting" << std::endl;
  os << *error_;
  Log::trace() << "ObsError<OBS>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSERROR_H_
