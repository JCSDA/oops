/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_LINEAROBSOPERATORS_H_
#define OOPS_BASE_LINEAROBSOPERATORS_H_

#include <string>

#include <boost/shared_ptr.hpp>

#include "util/Logger.h"
#include "oops/base/ObsSpace.h"
#include "oops/base/ObsOperators.h"
#include "oops/interface/LinearObsOperator.h"
#include "oops/interface/ModelAtLocations.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/interface/ObsVector.h"
#include "oops/interface/Variables.h"
#include "eckit/config/Configuration.h"
#include "util/DateTime.h"
#include "util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class LinearObsOperators : public util::Printable {
  typedef LinearObsOperator<MODEL>   LinearObsOperator_;
  typedef ModelAtLocations<MODEL>    ModelAtLocations_;
  typedef ObsAuxControl<MODEL>       ObsAuxControl_;
  typedef ObsAuxIncrement<MODEL>     ObsAuxIncrement_;
  typedef ObsOperators<MODEL>        ObsOperator_;
  typedef ObsVector<MODEL>           ObsVector_;
  typedef Variables<MODEL>           Variables_;

 public:
  static const std::string classname() {return "oops::LinearObsOperators";}

  explicit LinearObsOperators(const ObsOperator_ &);
  explicit LinearObsOperators(const LinearObsOperators &);
  ~LinearObsOperators();

/// Obs Operators
//  void setTrajectory(const ModelAtLocations_ &, const ObsAuxControl_ &);
//  void obsEquivTL(const ModelAtLocations_ &, ObsVector_ &, const ObsAuxIncrement_ &) const;
//  void obsEquivAD(ModelAtLocations_ &, const ObsVector_ &, ObsAuxIncrement_ &) const;

/// Other
  std::size_t size() const {return ops_.size();}
  LinearObsOperator_ & operator[](const std::size_t ii) {return *ops_.at(ii);} 
  const LinearObsOperator_ & operator[](const std::size_t ii) const {return *ops_.at(ii);} 
  Variables_ variables(const std::size_t) const;

 private:
  LinearObsOperators & operator=(const LinearObsOperators &);
  void print(std::ostream &) const;
  std::vector<boost::shared_ptr<LinearObsOperator_> > ops_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearObsOperators<MODEL>::LinearObsOperators(const ObsOperator_ & hop): ops_(0)
{
  for (std::size_t jobs = 0; jobs < hop.size(); ++jobs) {
    boost::shared_ptr<LinearObsOperator_> tmp(new LinearObsOperator_(hop[jobs]));
    ops_.push_back(tmp);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearObsOperators<MODEL>::LinearObsOperators(const LinearObsOperators & other)
 : ops_(other.ops_)
{}

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearObsOperators<MODEL>::~LinearObsOperators() {}

// -----------------------------------------------------------------------------

//template <typename MODEL>
//void LinearObsOperators<MODEL>::setTrajectory(const ModelAtLocations_ & gom,
//                                              const ObsAuxControl_ & aux) {
////  oper_->setTrajectory(gom.modelatlocations(), aux.obsauxcontrol());
//}

// -----------------------------------------------------------------------------

//template <typename MODEL>
//void LinearObsOperators<MODEL>::obsEquivTL(const ModelAtLocations_ & gom, ObsVector_ & yy,
//                                           const ObsAuxIncrement_ & aux) const {
////  oper_->obsEquivTL(gom.modelatlocations(), yy.obsvector(), aux.obsauxincrement());
//  for (std::size_t jobs = 0; jobs < ops_.size(); ++jobs) {
//    ops_[jobs]->obsEquivTL(gom[jobs], yy[jobs], aux);
//  }
//}

// -----------------------------------------------------------------------------

// template <typename MODEL>
// void LinearObsOperators<MODEL>::obsEquivAD(ModelAtLocations_ & gom, const ObsVector_ & yy,
//                                            ObsAuxIncrement_ & aux) const {
// //  oper_->obsEquivAD(gom.modelatlocations(), yy.obsvector(), aux.obsauxincrement());
//   for (std::size_t jobs = 0; jobs < ops_.size(); ++jobs) {
//     ops_[jobs]->obsEquivAD(gom[jobs], yy[jobs], aux);
//   }
// }

// -----------------------------------------------------------------------------

template <typename MODEL>
Variables<MODEL> LinearObsOperators<MODEL>::variables(const std::size_t jobs) const {
  Variables<MODEL> var(ops_.at(jobs)->variables());
  return var;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearObsOperators<MODEL>::print(std::ostream & os) const {
  for (std::size_t jobs = 0; jobs < ops_.size(); ++jobs) os << *ops_[jobs];
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_LINEAROBSOPERATORS_H_
