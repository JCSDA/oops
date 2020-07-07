/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#ifndef OOPS_BASE_QCDATA_H_
#define OOPS_BASE_QCDATA_H_

#include <vector>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include "oops/base/ObsSpaces.h"
#include "oops/interface/ObsDataVector.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief container for QC-related things (flags & obserrors)

template <typename OBS>
class QCData {
  typedef ObsSpaces<OBS>           ObsSpaces_;
  template <typename DATA> using ObsData_ = ObsDataVector<OBS, DATA>;
  template <typename DATA> using ObsDataPtr_ = boost::shared_ptr<ObsData_<DATA> >;

 public:
/// \brief Initializes QC data
  explicit QCData(const ObsSpaces_ &);

/// \brief accessor to QC flag
  const ObsDataPtr_<int> qcFlags(const size_t ii) const {return qcflags_[ii];}
/// \brief accessor to Obs errors
  const ObsDataPtr_<float> obsErrors(const size_t ii) const {return obserr_[ii];}

 private:
  std::vector<ObsDataPtr_<int> > qcflags_;   // QC flags
  std::vector<ObsDataPtr_<float> > obserr_;  // Obs Errors
};


// -----------------------------------------------------------------------------

template <typename OBS>
QCData<OBS>::QCData(const ObsSpaces_ & obspaces) {
  qcflags_.reserve(obspaces.size());
  obserr_.reserve(obspaces.size());
  for (size_t jj = 0; jj < obspaces.size(); ++jj) {
//  Allocate QC flags
    qcflags_.emplace_back(boost::make_shared<ObsData_<int>>(obspaces[jj],
                            obspaces[jj].obsvariables()));
//  Allocate and read initial obs error
    obserr_.emplace_back(boost::make_shared<ObsData_<float>>(obspaces[jj],
                            obspaces[jj].obsvariables(), "ObsError"));
  }
}

}  // namespace oops

#endif  // OOPS_BASE_QCDATA_H_
