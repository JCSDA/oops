/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_VARIABLECHANGEBASE_H_
#define OOPS_BASE_VARIABLECHANGEBASE_H_

#include <boost/noncopyable.hpp>

#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

/// Base class for generic variable transform

template <typename MODEL>
class VariableChangeBase : public util::Printable,
                           private boost::noncopyable {
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
  VariableChangeBase() {}
  virtual ~VariableChangeBase() {}

  virtual void linearize(const State_ &) =0;

  virtual void transform(const Increment_ &, Increment_ &) const =0;
  virtual void transformInverse(const Increment_ &, Increment_ &) const =0;
  virtual void transformAdjoint(const Increment_ &, Increment_ &) const =0;

 private:
  virtual void print(std::ostream &) const =0;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_VARIABLECHANGEBASE_H_
