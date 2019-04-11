/*
 * (C) Copyright 2018  UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_INTERFACE_OBSDATAVECTOR_H_
#define OOPS_INTERFACE_OBSDATAVECTOR_H_

#include <math.h>
#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/ObservationSpace.h"
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

  ObsDataVector(const ObservationSpace<MODEL> &, const Variables &);
  explicit ObsDataVector(const ObsDataVector &);
  ~ObsDataVector();

/// Interfacing
  ObsDataVec_ & obsdatavector() {return *data_;}
  const ObsDataVec_ & obsdatavector() const {return *data_;}

  ObsDataVector & operator = (const ObsDataVector &);

  void zero();
  unsigned int nobs() const {return data_->nobs();}

// I/O
  void read(const std::string &);
  void save(const std::string &) const;

 private:
  void print(std::ostream &) const;
  boost::scoped_ptr<ObsDataVec_> data_;
};

// -----------------------------------------------------------------------------
template <typename MODEL, typename DATATYPE>
ObsDataVector<MODEL, DATATYPE>::ObsDataVector(const ObservationSpace<MODEL> & os,
                                              const Variables & vars): data_() {
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::ObsDataVector starting" << std::endl;
  util::Timer timer(classname(), "ObsDataVector");
  data_.reset(new ObsDataVec_(os.observationspace(), vars));
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
void ObsDataVector<MODEL, DATATYPE>::print(std::ostream & os) const {
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *data_;
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::print done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename DATATYPE>
void ObsDataVector<MODEL, DATATYPE>::read(const std::string & name) {
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::read starting" << std::endl;
  util::Timer timer(classname(), "read");
  data_->read(name);
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::read done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename DATATYPE>
void ObsDataVector<MODEL, DATATYPE>::save(const std::string & name) const {
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::save starting";
  util::Timer timer(classname(), "save");
  data_->save(name);
  Log::trace() << "ObsDataVector<MODEL, DATATYPE>::save done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSDATAVECTOR_H_
