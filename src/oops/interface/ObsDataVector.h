/*
 * (C) Copyright 2018  UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_OBSDATAVECTOR_H_
#define OOPS_INTERFACE_OBSDATAVECTOR_H_

#include <math.h>
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

template <typename MODEL, typename DATATYPE>
class ObsDataVector : public util::Printable,
                      private util::ObjectCounter<ObsDataVector<MODEL, DATATYPE> > {
  typedef typename MODEL::template ObsDataVector<DATATYPE>  ObsDataVec_;

 public:
  static const std::string classname() {return "oops::ObsDataVector";}

  ObsDataVector(const ObsSpace<MODEL> &, const Variables &, const std::string name = "");
  explicit ObsDataVector(const ObsDataVector &);
  ~ObsDataVector();

/// Interfacing
  ObsDataVec_ & obsdatavector() {return *data_;}
  const ObsDataVec_ & obsdatavector() const {return *data_;}

  boost::shared_ptr<ObsDataVec_> obsdatavectorptr() {return data_;}
  boost::shared_ptr<const ObsDataVec_> obsdatavectorptr() const {return data_;}

  ObsDataVector & operator = (const ObsDataVector &);

  void zero();
  void mask(const ObsDataVector<MODEL, int> &);
  unsigned int nobs() const {return data_->nobs();}

// I/O
  void save(const std::string &) const;

 private:
  void print(std::ostream &) const;
  boost::shared_ptr<ObsDataVec_> data_;
};

// -----------------------------------------------------------------------------
template <typename MODEL, typename DATATYPE>
ObsDataVector<MODEL, DATATYPE>::ObsDataVector(const ObsSpace<MODEL> & os,
                                              const Variables & vars, const std::string name)
  : data_()
{
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::ObsDataVector starting" << std::endl;
  util::Timer timer(classname(), "ObsDataVector");
  data_.reset(new ObsDataVec_(os.obsspace(), vars, name));
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::ObsDataVector done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename DATATYPE>
ObsDataVector<MODEL, DATATYPE>::ObsDataVector(const ObsDataVector & other): data_() {
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::ObsDataVector starting" << std::endl;
  util::Timer timer(classname(), "ObsDataVector");
  data_.reset(new ObsDataVec_(*other.data_));
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::ObsDataVector done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename DATATYPE>
ObsDataVector<MODEL, DATATYPE>::~ObsDataVector() {
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::~ObsDataVector starting" << std::endl;
  util::Timer timer(classname(), "~ObsDataVector");
  data_.reset();
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::~ObsDataVector done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename DATATYPE> ObsDataVector<MODEL, DATATYPE> &
ObsDataVector<MODEL, DATATYPE>::operator=(const ObsDataVector & rhs) {
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  *data_ = *rhs.data_;
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::operator= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename DATATYPE>
void ObsDataVector<MODEL, DATATYPE>::zero() {
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");
  data_->zero();
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::zero done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename DATATYPE>
void ObsDataVector<MODEL, DATATYPE>::mask(const ObsDataVector<MODEL, int> & qc) {
  Log::trace() << "ObsDataVector<MODEL>::mask starting" << std::endl;
  util::Timer timer(classname(), "mask");
  data_->mask(qc.obsdatavector());
  Log::trace() << "ObsDataVector<MODEL>::mask done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename DATATYPE>
void ObsDataVector<MODEL, DATATYPE>::print(std::ostream & os) const {
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *data_;
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::print done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename DATATYPE>
void ObsDataVector<MODEL, DATATYPE>::save(const std::string & name) const {
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::save starting " << name << std::endl;
  util::Timer timer(classname(), "save");
  data_->save(name);
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::save done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
bool compareFlags(const ObsDataVector<MODEL, int> & first,
                  const ObsDataVector<MODEL, int> & second) {
  Log::trace() << "compareFlags(ObsDataVector<MODEL, int>) starting" << std::endl;
  bool compare = compareFlags(first.obsdatavector(), second.obsdatavector());
  Log::trace() << "compareFlags(ObsDataVector<MODEL, int>) done" << std::endl;
  return compare;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
size_t numZero(const ObsDataVector<MODEL, int> & data) {
  Log::trace() << "numZero(ObsDataVector<MODEL, int>) starting" << std::endl;
  size_t nzero = numZero(data.obsdatavector());
  Log::trace() << "numZero(ObsDataVector<MODEL, int>) done" << std::endl;
  return nzero;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSDATAVECTOR_H_
