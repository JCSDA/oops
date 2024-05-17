/*
 * (C) Copyright 2018  UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_OBSDATAVECTOR_H_
#define OOPS_INTERFACE_OBSDATAVECTOR_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/ObsVariables.h"
#include "oops/base/ObsVector.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "oops/interface/ObsDataVector_head.h"

namespace oops {

// -----------------------------------------------------------------------------
template <typename OBS, typename DATATYPE>
ObsDataVector<OBS, DATATYPE>::ObsDataVector(const ObsSpace<OBS> & os,
                                            const ObsVariables & vars, const std::string name)
  : data_()
{
  Log::trace() << "ObsDataVector<OBS, DATATYPE>::ObsDataVector starting" << std::endl;
  util::Timer timer(classname(), "ObsDataVector");
  data_.reset(new ObsDataVec_(os.obsspace(), vars, name));
  Log::trace() << "ObsDataVector<OBS, DATATYPE>::ObsDataVector done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS, typename DATATYPE>
ObsDataVector<OBS, DATATYPE>::ObsDataVector(const ObsDataVector & other): data_() {
  Log::trace() << "ObsDataVector<OBS, DATATYPE>::ObsDataVector starting" << std::endl;
  util::Timer timer(classname(), "ObsDataVector");
  data_.reset(new ObsDataVec_(*other.data_));
  Log::trace() << "ObsDataVector<OBS, DATATYPE>::ObsDataVector done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS, typename DATATYPE>
ObsDataVector<OBS, DATATYPE>::ObsDataVector(ObsVector<OBS> & other): data_() {
  Log::trace() << "ObsDataVector<OBS, DATATYPE>::ObsDataVector starting" << std::endl;
  util::Timer timer(classname(), "ObsDataVector");
  data_.reset(new ObsDataVec_(other.obsvector()));
  Log::trace() << "ObsDataVector<OBS, DATATYPE>::ObsDataVector done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS, typename DATATYPE>
ObsDataVector<OBS, DATATYPE>::~ObsDataVector() {
  Log::trace() << "ObsDataVector<OBS, DATATYPE>::~ObsDataVector starting" << std::endl;
  util::Timer timer(classname(), "~ObsDataVector");
  data_.reset();
  Log::trace() << "ObsDataVector<OBS, DATATYPE>::~ObsDataVector done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS, typename DATATYPE> ObsDataVector<OBS, DATATYPE> &
ObsDataVector<OBS, DATATYPE>::operator=(const ObsDataVector & rhs) {
  Log::trace() << "ObsDataVector<OBS, DATATYPE>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  *data_ = *rhs.data_;
  Log::trace() << "ObsDataVector<OBS, DATATYPE>::operator= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template <typename OBS, typename DATATYPE>
void ObsDataVector<OBS, DATATYPE>::zero() {
  Log::trace() << "ObsDataVector<OBS, DATATYPE>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");
  data_->zero();
  Log::trace() << "ObsDataVector<OBS, DATATYPE>::zero done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS, typename DATATYPE>
void ObsDataVector<OBS, DATATYPE>::mask(const ObsDataVector<OBS, int> & qc) {
  Log::trace() << "ObsDataVector<OBS>::mask starting" << std::endl;
  util::Timer timer(classname(), "mask");
  data_->mask(qc.obsdatavector());
  Log::trace() << "ObsDataVector<OBS>::mask done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS, typename DATATYPE>
void ObsDataVector<OBS, DATATYPE>::print(std::ostream & os) const {
  Log::trace() << "ObsDataVector<OBS, DATATYPE>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *data_;
  Log::trace() << "ObsDataVector<OBS, DATATYPE>::print done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS, typename DATATYPE>
void ObsDataVector<OBS, DATATYPE>::read(const std::string & name) {
  Log::trace() << "ObsDataVector<OBS, DATATYPE>::read starting " << name << std::endl;
  util::Timer timer(classname(), "read");
  data_->read(name);
  Log::trace() << "ObsDataVector<OBS, DATATYPE>::read done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS, typename DATATYPE>
void ObsDataVector<OBS, DATATYPE>::save(const std::string & name) const {
  Log::trace() << "ObsDataVector<OBS, DATATYPE>::save starting " << name << std::endl;
  util::Timer timer(classname(), "save");
  data_->save(name);
  Log::trace() << "ObsDataVector<OBS, DATATYPE>::save done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSDATAVECTOR_H_
