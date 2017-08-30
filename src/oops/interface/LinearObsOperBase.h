/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_LINEAROBSOPERBASE_H_
#define OOPS_INTERFACE_LINEAROBSOPERBASE_H_

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class LinearObsOperBase : public util::Printable,
                          private boost::noncopyable {
  typedef typename MODEL::ModelAtLocations    ModelAtLocations_;
  typedef typename MODEL::ObsAuxControl       ObsAuxControl_;
  typedef typename MODEL::ObsAuxIncrement     ObsAuxIncrement_;
  typedef typename MODEL::ObsVector           ObsVector_;
  typedef typename MODEL::Variables           Variables_;

 public:
  LinearObsOperBase() {}
  virtual ~LinearObsOperBase() {}

/// Obs Operators
  virtual void setTrajectory(const ModelAtLocations_ &, const ObsAuxControl_ &) =0;
  virtual void obsEquivTL(const ModelAtLocations_ &, ObsVector_ &, const ObsAuxIncrement_ &) const =0;
  virtual void obsEquivAD(ModelAtLocations_ &, const ObsVector_ &, ObsAuxIncrement_ &) const =0;

/// Other
  virtual boost::shared_ptr<const Variables_> variables() const =0;  // Required from Model

 private:
  virtual void print(std::ostream &) const =0;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_LINEAROBSOPERBASE_H_
