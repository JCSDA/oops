/*
 * (C) Copyright 2019  UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef LORENZ95_OBSDATA1D_H_
#define LORENZ95_OBSDATA1D_H_

#include <cmath>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/base/ObsVariables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "lorenz95/ObsTable.h"
#include "lorenz95/ObsVec1D.h"

namespace lorenz95 {

// -----------------------------------------------------------------------------
/// Data in observation space

template<typename DATATYPE>
class ObsData1D : public util::Printable,
                  private util::ObjectCounter<ObsData1D<DATATYPE> > {
 public:
  static const std::string classname() {return "lorenz95::ObsData1D";}

  ObsData1D(const ObsTable &, const oops::ObsVariables &, const std::string &);
  ObsData1D(const ObsData1D &);
  explicit ObsData1D(const ObsVec1D &);
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

  const ObsTable & obsdb_;
  std::vector<DATATYPE> data_;
};

//-----------------------------------------------------------------------------

template<typename DATATYPE>
ObsData1D<DATATYPE>::ObsData1D(const ObsTable & ot, const oops::ObsVariables &,
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
ObsData1D<DATATYPE>::ObsData1D(const ObsVec1D & other)
  : obsdb_(other.obsdb()), data_(other.size()) {
  const DATATYPE missing = util::missingValue<DATATYPE>();
  const double dmiss = util::missingValue<double>();
  for (size_t jj = 0; jj < data_.size(); ++jj) {
    if (other[jj] == dmiss) {
      data_.at(jj) = missing;
    } else {
      data_.at(jj) = static_cast<DATATYPE>(other[jj]);
    }
  }
}
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
  const DATATYPE missing = util::missingValue<DATATYPE>();
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
  const DATATYPE missing = util::missingValue<DATATYPE>();
  DATATYPE zmin = std::numeric_limits<DATATYPE>::max();
  DATATYPE zmax = std::numeric_limits<DATATYPE>::lowest();
  DATATYPE zavg = 0.0;
  size_t iobs = 0;
  for (const DATATYPE & val : data_) {
    if (val != missing) {
      if (val < zmin) zmin = val;
      if (val > zmax) zmax = val;
      zavg += val;
      ++iobs;
    }
  }
  if (iobs > 0) {
    zavg /= static_cast<DATATYPE>(iobs);
    os << "Lorenz 95 nobs= " << iobs << " Min=" << zmin << ", Max=" << zmax
       << ", Average=" << zavg;
  } else {
    os << "Lorenz 95 : No observations";
  }
}
// -----------------------------------------------------------------------------
}  // namespace lorenz95

#endif  // LORENZ95_OBSDATA1D_H_
