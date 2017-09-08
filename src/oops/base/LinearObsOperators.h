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

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "oops/base/ObsSpaces.h"
#include "oops/interface/LinearObsOperator.h"
#include "oops/interface/ModelAtLocations.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/interface/ObsVector.h"
#include "oops/interface/Variables.h"
#include "util/DateTime.h"
#include "util/Logger.h"
#include "util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class LinearObsOperators : public util::Printable,
                           private boost::noncopyable {
  typedef LinearObsOperator<MODEL>   LinearObsOperator_;
  typedef ModelAtLocations<MODEL>    ModelAtLocations_;
  typedef ObsAuxControl<MODEL>       ObsAuxControl_;
  typedef ObsAuxIncrement<MODEL>     ObsAuxIncrement_;
  typedef ObsSpaces<MODEL>           ObsSpace_;
  typedef ObsVector<MODEL>           ObsVector_;
  typedef Variables<MODEL>           Variables_;

 public:
  static const std::string classname() {return "oops::LinearObsOperators";}

  explicit LinearObsOperators(const ObsSpace_ &);
  ~LinearObsOperators();

/// Access
  std::size_t size() const {return ops_.size();}
  LinearObsOperator_ & operator[](const std::size_t ii) {return *ops_.at(ii);}
  const LinearObsOperator_ & operator[](const std::size_t ii) const {return *ops_.at(ii);}
  Variables_ variables(const std::size_t) const;

 private:
  void print(std::ostream &) const;
  std::vector<boost::shared_ptr<LinearObsOperator_> > ops_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearObsOperators<MODEL>::LinearObsOperators(const ObsSpace_ & os): ops_(0) {
  for (std::size_t jobs = 0; jobs < os.size(); ++jobs) {
    boost::shared_ptr<LinearObsOperator_> tmp(new LinearObsOperator_(os[jobs]));
    ops_.push_back(tmp);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearObsOperators<MODEL>::~LinearObsOperators() {}

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
