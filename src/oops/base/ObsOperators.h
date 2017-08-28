/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_OBSOPERATORS_H_
#define OOPS_BASE_OBSOPERATORS_H_

#include <string>

#include <boost/shared_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/ObsSpace.h"
#include "oops/interface/ModelAtLocations.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsVector.h"
#include "oops/interface/Variables.h"
#include "util/DateTime.h"
#include "util/Logger.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsOperators : public util::Printable,
                     private util::ObjectCounter<ObsOperators<MODEL> > {
  typedef ModelAtLocations<MODEL>    ModelAtLocations_;
  typedef ObsAuxControl<MODEL>       ObsAuxControl_;
  typedef ObsOperator<MODEL>         ObsOperator_;
  typedef ObsVector<MODEL>           ObsVector_;
  typedef ObsSpace<MODEL>            ObsSpace_;
  typedef Variables<MODEL>           Variables_;

 public:
  static const std::string classname() {return "oops::ObsOperators";}

  explicit ObsOperators(const ObsSpace_ &);
  ObsOperators(const ObsOperators &);
  ~ObsOperators();

/// Obs Operator
//  void obsEquiv(const ModelAtLocations_ &, ObsVector_ &, const ObsAuxControl_ &) const;

/// Other
  std::size_t size() const {return ops_.size();}
  const ObsOperator_ & operator[](const std::size_t ii) const {return *ops_.at(ii);} 
  Variables_ variables(const std::size_t jobs) const;

 private:
  ObsOperators & operator=(const ObsOperators &);
  void print(std::ostream &) const;
  std::vector<boost::shared_ptr<ObsOperator_> > ops_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsOperators<MODEL>::ObsOperators(const ObsSpace_ & os) : ops_(0)
{
  for (std::size_t jobs = 0; jobs < os.size(); ++jobs) {
    boost::shared_ptr<ObsOperator_> tmp(new ObsOperator_(os[jobs]));
    ops_.push_back(tmp);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsOperators<MODEL>::ObsOperators(const ObsOperators & other) : ops_(other.ops_)
{}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsOperators<MODEL>::~ObsOperators() {}

// -----------------------------------------------------------------------------

// template <typename MODEL>
// void ObsOperators<MODEL>::obsEquiv(const ModelAtLocations_ & gom, ObsVector_ & yy,
//                                    const ObsAuxControl_ & aux) const {
//   for (std::size_t jobs = 0; jobs < ops_.size(); ++jobs) {
//     ops_[jobs]->obsEquiv(gom[jobs], yy[jobs], aux);
//   }
// }

// -----------------------------------------------------------------------------

template <typename MODEL>
Variables<MODEL> ObsOperators<MODEL>::variables(const std::size_t jobs) const {
  Variables<MODEL> var(ops_.at(jobs)->variables());
  return var;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsOperators<MODEL>::print(std::ostream & os) const {
  for (std::size_t jobs = 0; jobs < ops_.size(); ++jobs) os << *ops_[jobs];
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSOPERATORS_H_
