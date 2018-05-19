/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_OBSERRORS_H_
#define OOPS_BASE_OBSERRORS_H_

#include <string>
#include <vector>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Departures.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/interface/ObsErrorCovariance.h"
#include "oops/interface/ObsVector.h"
#include "util/Logger.h"
#include "util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsErrors : public util::Printable,
                  private boost::noncopyable {
  typedef Departures<MODEL>          Departures_;
  typedef Observations<MODEL>        Observations_;
  typedef ObsErrorCovariance<MODEL>  ObsError_;
  typedef ObsSpaces<MODEL>           ObsSpace_;
  typedef ObsVector<MODEL>           ObsVector_;

 public:
  static const std::string classname() {return "oops::ObsErrors";}

  explicit ObsErrors(const ObsSpace_ &);
  ~ObsErrors();

/// Access
  std::size_t size() const {return err_.size();}
  const ObsError_ & operator[](const std::size_t ii) const {return *err_.at(ii);}

/// Linearize and reset for inner loop if needed
  void linearize(const Observations_ &);

/// Multiply a Departure by \f$R\f$ and \f$R^{-1}\f$
  Departures_ * multiply(const Departures_ &) const;
  Departures_ * inverseMultiply(const Departures_ &) const;

/// Generate random perturbation
  void randomize(Departures_ &) const;

 private:
  void print(std::ostream &) const;
  std::vector<boost::shared_ptr<ObsError_> > err_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsErrors<MODEL>::ObsErrors(const ObsSpace_ & os) : err_(0)
{
  for (std::size_t jj = 0; jj < os.size(); ++jj) {
    eckit::LocalConfiguration conf(os[jj].config(), "Covariance");
    boost::shared_ptr<ObsError_> tmp(new ObsError_(os[jj], conf));
    err_.push_back(tmp);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsErrors<MODEL>::~ObsErrors() {}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObsErrors<MODEL>::linearize(const Observations_ & yy) {
  for (std::size_t jj = 0; jj < err_.size(); ++jj) {
    err_[jj]->linearize(yy[jj]);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Departures<MODEL> * ObsErrors<MODEL>::multiply(const Departures_ & dy) const {
  std::vector<boost::shared_ptr<ObsVector_> > ovec;
  for (std::size_t jj = 0; jj < err_.size(); ++jj) {
    boost::shared_ptr<ObsVector_> tmp(err_[jj]->multiply(dy[jj]));
    ovec.push_back(tmp);
  }
  return new Departures_(ovec);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Departures<MODEL> * ObsErrors<MODEL>::inverseMultiply(const Departures_ & dy) const {
  std::vector<boost::shared_ptr<ObsVector_> > ovec;
  for (std::size_t jj = 0; jj < err_.size(); ++jj) {
    boost::shared_ptr<ObsVector_> tmp(err_[jj]->inverseMultiply(dy[jj]));
    ovec.push_back(tmp);
  }
  return new Departures_(ovec);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObsErrors<MODEL>::randomize(Departures_ & dy) const {
  for (std::size_t jj = 0; jj < err_.size(); ++jj) {
    err_[jj]->randomize(dy[jj]);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsErrors<MODEL>::print(std::ostream & os) const {
  for (std::size_t jj = 0; jj < err_.size(); ++jj) os << *err_[jj];
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERRORS_H_
