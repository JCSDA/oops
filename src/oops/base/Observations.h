/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_OBSERVATIONS_H_
#define OOPS_BASE_OBSERVATIONS_H_

#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "util/Logger.h"
#include "oops/base/Departures.h"
#include "oops/base/ObsSpace.h"
#include "oops/interface/ModelAtLocations.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsVector.h"
#include "util/DateTime.h"
#include "util/Printable.h"
#include "util/abor1_cpp.h"

namespace oops {

template<typename MODEL> class Observations;

/// Observations Class.
/*!
 *  Contains observed values or their model equivalents and the associated
 *  observation operator.
 */

// -----------------------------------------------------------------------------
template <typename MODEL> class Observations : public util::Printable {
  typedef Departures<MODEL>          Departures_;
  typedef ModelAtLocations<MODEL>    GOM_;
  typedef ObsAuxControl<MODEL>       ObsAuxCtrl_;
  typedef ObsOperator<MODEL>         ObsOperator_;
  typedef ObsSpace<MODEL>            ObsSpace_;
  typedef ObsVector<MODEL>           ObsVector_;

 public:
  explicit Observations(const ObsSpace_ &);
  explicit Observations(const Observations &);
  ~Observations();
  Observations & operator=(const Observations &);

/// Access
  std::size_t size() const {return obs_.size();}
  ObsVector_ & operator[](const std::size_t ii) {return obs_.at(ii);} 
  const ObsVector_ & operator[](const std::size_t ii) const {return obs_.at(ii);} 

/// Interactions with Departures
  std::vector<boost::shared_ptr<ObsVector_> > operator-(const Observations & other) const;
  Observations & operator+=(const Departures_ &);

/// Get observations values
  const ObsVector_ & obsvalues() const {return obs_;}

/// Save observations values
  void save(const std::string &) const;
  void read(const eckit::Configuration &);

 private:
  void print(std::ostream &) const;

  boost::ptr_vector<ObsVector_> obs_;
};

// =============================================================================

template <typename MODEL>
Observations<MODEL>::Observations(const ObsSpace_ & obsdb): obs_(0)
{
  for (std::size_t jj = 0; jj < obsdb.size(); ++jj) {
    obs_.push_back(new ObsVector_(obsdb[jj]));
  }
  Log::trace() << "Observations created" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Observations<MODEL>::Observations(const Observations & other): obs_(other.obs_)
{
  Log::trace() << "Observations copy-created" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Observations<MODEL>::~Observations() {
  Log::trace() << "Observations destructed" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Observations<MODEL> & Observations<MODEL>::operator=(const Observations & rhs) {
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) {
    obs_[jj] = rhs.obs_[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
std::vector<boost::shared_ptr<ObsVector<MODEL> > >
Observations<MODEL>::operator-(const Observations & other) const {
  std::vector<boost::shared_ptr<ObsVector_> > out;
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) {
    boost::shared_ptr<ObsVector_> ovec(new ObsVector_(obs_[jj], true));
    *ovec -= other.obs_[jj];
    out.push_back(ovec);
  }
  return out;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Observations<MODEL> & Observations<MODEL>::operator+=(const Departures_ & dy) {
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) {
    obs_[jj] += dy[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Observations<MODEL>::save(const std::string & name) const {
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) {
    obs_[jj].save(name);
  }
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Observations<MODEL>::read(const eckit::Configuration & config) {
  std::vector<eckit::LocalConfiguration> conf;
  config.get("ObsTypes", conf);
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) {
    const std::string name = conf[jj].getString("ObsData.obsvalue");
    obs_[jj].read(name);
  }
  Log::trace() << "Observations:Observations have been read" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Observations<MODEL>::print(std::ostream & os) const {
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) os << obs_[jj];
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_OBSERVATIONS_H_
