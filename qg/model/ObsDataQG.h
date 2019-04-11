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

#include "model/QgFortran.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace qg {
  class ObsSpaceQG;

// -----------------------------------------------------------------------------
/// Data in observation space

template<typename DATATYPE>
class ObsDataQG : public util::Printable,
                  private util::ObjectCounter<ObsDataQG<DATATYPE> > {
 public:
  static const std::string classname() {return "qg::ObsDataQG";}

  ObsDataQG(const ObsSpaceQG &, const oops::Variables &);
  ObsDataQG(const ObsDataQG &);
  ~ObsDataQG() {}

  ObsDataQG & operator= (const ObsDataQG &);

  void zero();

// I/O
  void read(const std::string &);
  void save(const std::string &) const;

 private:
  void print(std::ostream &) const;

//  const ObsSpaceQG & obsdb_;
};
//-----------------------------------------------------------------------------

template<typename DATATYPE>
ObsDataQG<DATATYPE>::ObsDataQG(const ObsSpaceQG &, const oops::Variables &) {
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
ObsDataQG<DATATYPE>::ObsDataQG(const ObsDataQG &) {
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
ObsDataQG<DATATYPE> & ObsDataQG<DATATYPE>::operator= (const ObsDataQG &) {
  return *this;
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsDataQG<DATATYPE>::zero() {
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsDataQG<DATATYPE>::read(const std::string &) {
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsDataQG<DATATYPE>::save(const std::string &) const {
}
// -----------------------------------------------------------------------------
template<typename DATATYPE>
void ObsDataQG<DATATYPE>::print(std::ostream &) const {
}
// -----------------------------------------------------------------------------
}  // namespace qg

#endif  // QG_MODEL_OBSDATAQG_H_
