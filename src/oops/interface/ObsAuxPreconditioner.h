/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_OBSAUXPRECONDITIONER_H_
#define OOPS_INTERFACE_OBSAUXPRECONDITIONER_H_

#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include "oops/interface/ObsAuxIncrement.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"


namespace oops {

// -----------------------------------------------------------------------------
/// \brief This class calls the <OBS> implementation of ObsAuxPreconditioner.
// -----------------------------------------------------------------------------

template <typename OBS>
class ObsAuxPreconditioner : public util::Printable,
                         private util::ObjectCounter<ObsAuxPreconditioner<OBS> > {
  typedef typename OBS::ObsAuxPreconditioner      ObsAuxPreconditioner_;
  typedef ObsAuxIncrement<OBS>                    ObsAuxIncrement_;

 public:
  static const std::string classname() {return "oops::ObsAuxPreconditioner";}

  /// Constructor
  explicit ObsAuxPreconditioner(std::unique_ptr<ObsAuxPreconditioner_>);
  /// Destructor
  ~ObsAuxPreconditioner();
  ObsAuxPreconditioner(ObsAuxPreconditioner &&) = default;
  ObsAuxPreconditioner & operator = (ObsAuxPreconditioner &&) = default;

  /// Sets the second parameter to the first multiplied by the preconditioner matrix.
  void multiply(const ObsAuxIncrement_ &, ObsAuxIncrement_ &) const;

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ObsAuxPreconditioner_> precon_;
};

// =============================================================================

template<typename OBS>
ObsAuxPreconditioner<OBS>::ObsAuxPreconditioner(std::unique_ptr<ObsAuxPreconditioner_> precon)
    : precon_(std::move(precon))
{
  Log::trace() << "ObsAuxPreconditioner<OBS>::ObsAuxPreconditioner constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
ObsAuxPreconditioner<OBS>::~ObsAuxPreconditioner() {
  Log::trace() << "ObsAuxPreconditioner<OBS>::~ObsAuxPreconditioner starting" << std::endl;
  util::Timer timer(classname(), "~ObsAuxPreconditioner");
  precon_.reset();
  Log::trace() << "ObsAuxPreconditioner<OBS>::~ObsAuxPreconditioner done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsAuxPreconditioner<OBS>::multiply(const ObsAuxIncrement_ & dx1,
                                         ObsAuxIncrement_ & dx2) const {
  Log::trace() << "ObsAuxPreconditioner<OBS>::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  precon_->multiply(dx1.obsauxincrement(), dx2.obsauxincrement());
  Log::trace() << "ObsAuxPreconditioner<OBS>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsAuxPreconditioner<OBS>::print(std::ostream & os) const {
  Log::trace() << "ObsAuxPreconditioner<OBS>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *precon_;
  Log::trace() << "ObsAuxPreconditioner<OBS>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSAUXPRECONDITIONER_H_
