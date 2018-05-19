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
#include <vector>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "oops/base/ObsSpaces.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsOperators : public util::Printable,
                     private boost::noncopyable {
  typedef ObsAuxControl<MODEL>       ObsAuxControl_;
  typedef ObsOperator<MODEL>         ObsOperator_;
  typedef ObsVector<MODEL>           ObsVector_;
  typedef ObsSpaces<MODEL>           ObsSpace_;

 public:
  static const std::string classname() {return "oops::ObsOperators";}

  explicit ObsOperators(const ObsSpace_ &);
  ~ObsOperators();

/// Access
  std::size_t size() const {return ops_.size();}
  const ObsOperator_ & operator[](const std::size_t ii) const {return *ops_.at(ii);}
  const Variables & variables(const std::size_t jobs) const;

 private:
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
ObsOperators<MODEL>::~ObsOperators() {}

// -----------------------------------------------------------------------------

template <typename MODEL>
const Variables & ObsOperators<MODEL>::variables(const std::size_t jobs) const {
  return ops_.at(jobs)->variables();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsOperators<MODEL>::print(std::ostream & os) const {
  for (std::size_t jobs = 0; jobs < ops_.size(); ++jobs) os << *ops_[jobs];
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSOPERATORS_H_
