/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSAUXPRECONDITIONERS_H_
#define OOPS_BASE_OBSAUXPRECONDITIONERS_H_

#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "oops/base/ObsAuxIncrements.h"
#include "oops/interface/ObsAuxPreconditioner.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Holds a vector of ObsAuxPreconditioner
// -----------------------------------------------------------------------------

template <typename OBS>
class ObsAuxPreconditioners : public util::Printable {
  typedef ObsAuxPreconditioner<OBS>    ObsAuxPreconditioner_;
  typedef ObsAuxIncrements<OBS>        ObsAuxIncrements_;
  typedef ObsSpaces<OBS>               ObsSpaces_;

 public:
  static const std::string classname() {return "oops::ObsAuxPreconditioners";}

  explicit ObsAuxPreconditioners(std::vector<ObsAuxPreconditioner<OBS>>);
  ~ObsAuxPreconditioners() = default;
  ObsAuxPreconditioners(ObsAuxPreconditioners &&) = default;
  ObsAuxPreconditioners & operator = (ObsAuxPreconditioners &&) = default;

  /// Operators
  void multiply(const ObsAuxIncrements_ &, ObsAuxIncrements_ &) const;

 private:
  void print(std::ostream &) const;
  std::vector<ObsAuxPreconditioner_> precond_;
};

// =============================================================================

template<typename OBS>
ObsAuxPreconditioners<OBS>::ObsAuxPreconditioners(std::vector<ObsAuxPreconditioner<OBS>> precond)
  : precond_(std::move(precond))
{
  Log::trace() << "ObsAuxPreconditioners<OBS>::ObsAuxPreconditioners done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsAuxPreconditioners<OBS>::multiply(const ObsAuxIncrements_ & dx1,
                                       ObsAuxIncrements_ & dx2) const {
  Log::trace() << "ObsAuxPreconditioners<OBS>::multiply starting" << std::endl;
  ASSERT(dx1.size() == dx2.size() && precond_.size() == dx1.size());
  for (std::size_t jobs = 0; jobs < dx1.size(); ++jobs) {
    precond_[jobs].multiply(dx1[jobs], dx2[jobs]);
  }
  Log::trace() << "ObsAuxPreconditioners<OBS>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsAuxPreconditioners<OBS>::print(std::ostream & os) const {
  Log::trace() << "ObsAuxPreconditioners<OBS>::print starting" << std::endl;
  for (std::size_t jobs = 0; jobs < precond_.size(); ++jobs) os << precond_.at(jobs) << " ";
  Log::trace() << "ObsAuxPreconditioners<OBS>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSAUXPRECONDITIONERS_H_
