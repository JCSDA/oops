/*
 * (C) Copyright 2019  UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef LORENZ95_OBSDATA1D_H_
#define LORENZ95_OBSDATA1D_H_

#include <cmath>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "lorenz95/ObsTableView.h"

namespace lorenz95 {

// -----------------------------------------------------------------------------
/// Data in observation space

template<typename DATATYPE>
class ObsData1D : public util::Printable,
                  private util::ObjectCounter<ObsData1D<DATATYPE> > {
 public:
  static const std::string classname() {return "lorenz95::ObsData1D";}

  ObsData1D(const ObsTableView &, const oops::Variables &, const std::string &);
  ObsData1D(const ObsData1D &);
  ~ObsData1D() {}

  ObsData1D & operator= (const ObsData1D &);

  void zero();
  void mask(const ObsData1D<int> &);

  size_t nobs() const {return data_.size();}
  DATATYPE & operator[] (const size_t ii) {return data_.at(ii);}
  const DATATYPE & operator[] (const size_t ii) const {return data_.at(ii);}

// I/O
  void read(const std::string &);
  void save(const std::string &) const;

 private:
  void print(std::ostream &) const;

  const ObsTableView & obsdb_;
  std::vector<DATATYPE> data_;
};

//-----------------------------------------------------------------------------

template<typename DATATYPE>
ObsData1D<DATATYPE>::ObsData1D(const ObsTableView & ot, const oops::Variables &,
                               const std::string & name)
  : obsdb_(ot), data_(ot.nobs())
{
  this->zero();
  if (!name.empty()) obsdb_.getdb(name, data_);
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
ObsData1D<DATATYPE>::ObsData1D(const ObsData1D & other)
  : obsdb_(other.obsdb_), data_(other.data_)
{}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
ObsData1D<DATATYPE> & ObsData1D<DATATYPE>::operator= (const ObsData1D & rhs) {
  ASSERT(data_.size() == rhs.data_.size());
  data_ = rhs.data_;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsData1D<DATATYPE>::zero() {
  for (size_t jj = 0; jj < data_.size(); ++jj) {
    data_.at(jj) = static_cast<DATATYPE>(0);
  }
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsData1D<DATATYPE>::mask(const ObsData1D<int> & mask) {
  DATATYPE missing = util::missingValue(missing);
  for (size_t jj = 0; jj < data_.size(); ++jj) {
    if (mask[jj]) data_.at(jj) = missing;
  }
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsData1D<DATATYPE>::read(const std::string & name) {
  obsdb_.getdb(name, data_);
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsData1D<DATATYPE>::save(const std::string & name) const {
  obsdb_.putdb(name, data_);
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsData1D<DATATYPE>::print(std::ostream & os) const {
  if (data_.size() > 0) {
    DATATYPE zmin = data_.at(0);
    DATATYPE zmax = data_.at(0);
    for (size_t jj = 0; jj < data_.size(); ++jj) {
      if (data_.at(jj) < zmin) zmin = data_.at(jj);
      if (data_.at(jj) > zmax) zmax = data_.at(jj);
    }
    os << "Lorenz 95 nobs= " << data_.size() << " Min=" << zmin << ", Max=" << zmax;
  } else {
    os << "Lorenz 95 nobs= " << data_.size() << " --- No observations";
  }
}
// -----------------------------------------------------------------------------
}  // namespace lorenz95

#endif  // LORENZ95_OBSDATA1D_H_
