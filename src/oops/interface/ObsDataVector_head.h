/*
 * (C) Copyright 2018  UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_OBSDATAVECTOR_HEAD_H_
#define OOPS_INTERFACE_OBSDATAVECTOR_HEAD_H_

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
  template <typename OBS> class ObsVector;

// -----------------------------------------------------------------------------

template <typename OBS, typename DATATYPE>
class ObsDataVector : public util::Printable,
                      private util::ObjectCounter<ObsDataVector<OBS, DATATYPE> > {
  typedef typename OBS::template ObsDataVector<DATATYPE>  ObsDataVec_;

 public:
  static const std::string classname() {return "oops::ObsDataVector";}

  ObsDataVector(const ObsSpace<OBS> &, const Variables &, const std::string name = "");
  explicit ObsDataVector(const ObsDataVector &);
  explicit ObsDataVector(ObsVector<OBS> &);
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

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSDATAVECTOR_HEAD_H_
