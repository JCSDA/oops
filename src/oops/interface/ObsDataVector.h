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

#include "oops/base/Variables.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename OBS, typename DATATYPE>
class ObsDataVector : public util::Printable,
                      private util::ObjectCounter<ObsDataVector<OBS, DATATYPE> > {
  typedef typename OBS::template ObsDataVector<DATATYPE>  ObsDataVec_;

 public:
  static const std::string classname() {return "oops::ObsDataVector";}

  ObsDataVector(const ObsSpace<OBS> &, const Variables &, const std::string name = "");
  explicit ObsDataVector(const ObsDataVector &);
  ~ObsDataVector();

/// Interfacing
  ObsDataVec_ & obsdatavector() {return *data_;}
  const ObsDataVec_ & obsdatavector() const {return *data_;}

  std::shared_ptr<ObsDataVec_> obsdatavectorptr() {return data_;}
  std::shared_ptr<const ObsDataVec_> obsdatavectorptr() const {return data_;}

  ObsDataVector & operator = (const ObsDataVector &);

  void zero();
  void mask(const ObsDataVector<OBS, int> &);
  unsigned int nobs() const {return data_->nobs();}

// I/O
  void read(const std::string &);
  void save(const std::string &) const;

 private:
  void print(std::ostream &) const;
  std::shared_ptr<ObsDataVec_> data_;
};

// -----------------------------------------------------------------------------
template <typename OBS, typename DATATYPE>
ObsDataVector<OBS, DATATYPE>::ObsDataVector(const ObsSpace<OBS> & os,
                                              const Variables & vars, const std::string name)
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
