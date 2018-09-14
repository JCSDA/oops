/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSENSEMBLE_H_
#define OOPS_BASE_OBSENSEMBLE_H_

#include <boost/ptr_container/ptr_vector.hpp>

#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Ensemble

template<typename MODEL> class ObsEnsemble {
  typedef Observations<MODEL>        Observations_;
  typedef ObsSpaces<MODEL>           ObsSpace_;

 public:
/// Constructor
  ObsEnsemble(const ObsSpace_ &, const int &);

/// Destructor
  virtual ~ObsEnsemble() {}

  /// Accessors
  unsigned int size() const {
    return rank_;
  }
  Observations_ & operator[](const int ii) {
    return ensemblePerturbs_[ii];
  }
  const Observations_ & operator[](const int ii) const {
    return ensemblePerturbs_[ii];
  }

 private:
  unsigned int rank_;
  boost::ptr_vector<Observations_> ensemblePerturbs_;
};

// ====================================================================================

template<typename MODEL>
ObsEnsemble<MODEL>::ObsEnsemble(const ObsSpace_ & os, const int & rank)
  : rank_(rank),
    ensemblePerturbs_()
{
  for (unsigned i = 0; i < rank_; ++i) {
    Observations_ * y = new Observations_(os);
    ensemblePerturbs_.push_back(y);
  }
  Log::trace() << "ObsEnsemble:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSENSEMBLE_H_
