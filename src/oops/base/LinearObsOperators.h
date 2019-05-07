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
#include <vector>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "oops/base/ObsSpaces.h"
#include "oops/base/Variables.h"
#include "oops/interface/LinearObsOperator.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class LinearObsOperators : public util::Printable,
                           private boost::noncopyable {
  typedef LinearObsOperator<MODEL>   LinearObsOperator_;
  typedef ObsAuxControl<MODEL>       ObsAuxControl_;
  typedef ObsAuxIncrement<MODEL>     ObsAuxIncrement_;
  typedef ObsSpaces<MODEL>           ObsSpace_;
  typedef ObsVector<MODEL>           ObsVector_;

 public:
  static const std::string classname() {return "oops::LinearObsOperators";}

  LinearObsOperators(const ObsSpace_ &, const eckit::Configuration &);
  ~LinearObsOperators();

/// Access
  std::size_t size() const {return ops_.size();}
  LinearObsOperator_ & operator[](const std::size_t ii) {return *ops_.at(ii);}
  const LinearObsOperator_ & operator[](const std::size_t ii) const {return *ops_.at(ii);}
  const Variables & variables(const std::size_t) const;

 private:
  void print(std::ostream &) const;
  std::vector<boost::shared_ptr<LinearObsOperator_> > ops_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearObsOperators<MODEL>::LinearObsOperators(const ObsSpace_ & os,
                                              const eckit::Configuration & conf) : ops_(0) {
  std::vector<eckit::LocalConfiguration> typeconfs;
  conf.get("ObsTypes", typeconfs);
  for (std::size_t jobs = 0; jobs < os.size(); ++jobs) {
    eckit::LocalConfiguration obsopconf(typeconfs[jobs], "ObsOperator");
    boost::shared_ptr<LinearObsOperator_> tmp(new LinearObsOperator_(os[jobs], obsopconf));
    ops_.push_back(tmp);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearObsOperators<MODEL>::~LinearObsOperators() {}

// -----------------------------------------------------------------------------

template <typename MODEL>
const Variables & LinearObsOperators<MODEL>::variables(const std::size_t jobs) const {
  return ops_.at(jobs)->variables();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearObsOperators<MODEL>::print(std::ostream & os) const {
  for (std::size_t jobs = 0; jobs < ops_.size(); ++jobs) os << *ops_[jobs];
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_LINEAROBSOPERATORS_H_
