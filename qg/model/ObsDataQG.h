/*
 * (C) Copyright 2019  UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef QG_MODEL_OBSDATAQG_H_
#define QG_MODEL_OBSDATAQG_H_

#include <math.h>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/ObsVariables.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "oops/qg/ObsSpaceQG.h"
#include "oops/qg/ObsVecQG.h"
#include "oops/qg/QgFortran.h"

namespace qg {

// -----------------------------------------------------------------------------
/// Data in observation space

template<typename DATATYPE>
class ObsDataQG : public util::Printable,
                  private util::ObjectCounter<ObsDataQG<DATATYPE> > {
 public:
  static const std::string classname() {return "qg::ObsDataQG";}

  explicit ObsDataQG(const ObsSpaceQG &, const oops::ObsVariables &,
                     const std::string &);
  ObsDataQG(const ObsDataQG &);
  explicit ObsDataQG(const ObsVecQG &);
  ~ObsDataQG() {}

  ObsDataQG & operator= (const ObsDataQG &);

  /// set all values to zero
  void zero();
  /// set all values to one
  void ones();
  void mask(const ObsDataQG<int>);

// I/O
  void read(const std::string &);
  void save(const std::string &) const;

  const int & toFortran() const {return data_.toFortran();}
  const ObsVecQG & vect() const {return data_;}
 private:
  void print(std::ostream &) const;
  ObsVecQG data_;
};
//-----------------------------------------------------------------------------

template<typename DATATYPE>
ObsDataQG<DATATYPE>::ObsDataQG(const ObsSpaceQG & os, const oops::ObsVariables & var,
                               const std::string & name): data_(os, name) {
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
ObsDataQG<DATATYPE>::ObsDataQG(const ObsDataQG & other): data_(other.data_) {
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
ObsDataQG<DATATYPE>::ObsDataQG(const ObsVecQG & other): data_(other) {
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
ObsDataQG<DATATYPE> & ObsDataQG<DATATYPE>::operator= (const ObsDataQG & rhs) {
  data_ = rhs.data_;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsDataQG<DATATYPE>::zero() {
  data_.zero();
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsDataQG<DATATYPE>::ones() {
  data_.ones();
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsDataQG<DATATYPE>::mask(const ObsDataQG<int> mask) {
  qg_obsvec_mask_f90(data_.toFortran(), mask.toFortran());
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsDataQG<DATATYPE>::read(const std::string & name) {
  data_.read(name);
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsDataQG<DATATYPE>::save(const std::string & name) const {
  data_.save(name);
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsDataQG<DATATYPE>::print(std::ostream & os) const {
  os << data_;
}
// -----------------------------------------------------------------------------
}  // namespace qg

#endif  // QG_MODEL_OBSDATAQG_H_
