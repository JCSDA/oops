/*
 * (C) Copyright 2017-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSAUXCONTROLS_H_
#define OOPS_BASE_OBSAUXCONTROLS_H_

#include <iostream>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "oops/base/ObsSpaces.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsAuxControls : public util::Printable {
  typedef ObsAuxControl<MODEL>       ObsAuxControl_;
  typedef ObsSpaces<MODEL>           ObsSpaces_;

 public:
  static const std::string classname() {return "oops::ObsAuxControls";}

  ObsAuxControls(const ObsSpaces_ &, const eckit::Configuration &);
  explicit ObsAuxControls(const ObsAuxControls &, const bool copy = true);
  ~ObsAuxControls();

/// Access
  std::size_t size() const {return auxs_.size();}
  const ObsAuxControl_ & operator[](const std::size_t ii) const {return *auxs_.at(ii);}
  ObsAuxControl_ & operator[](const std::size_t ii) {return *auxs_.at(ii);}

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;

  ObsAuxControls & operator=(const ObsAuxControls &);

 private:
  void print(std::ostream &) const;
  std::vector<boost::shared_ptr<ObsAuxControl_> > auxs_;
};

// =============================================================================

template<typename MODEL>
ObsAuxControls<MODEL>::ObsAuxControls(const ObsSpaces_ & odb, const eckit::Configuration & conf)
  : auxs_(0)
{
  std::vector<eckit::LocalConfiguration> obsconf;
  conf.get("ObsTypes", obsconf);
  for (std::size_t jobs = 0; jobs < obsconf.size(); ++jobs) {
    boost::shared_ptr<ObsAuxControl_> tmp(new ObsAuxControl_(odb[jobs], obsconf[jobs]));
    auxs_.push_back(tmp);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ObsAuxControls<MODEL>::ObsAuxControls(const ObsAuxControls & other, const bool copy)
  : auxs_(other.size())
{
  Log::trace() << "ObsAuxControls<MODEL>::ObsAuxControls copy starting" << std::endl;
  for (std::size_t jobs = 0; jobs < other.size(); ++jobs) {
    auxs_[jobs].reset(new ObsAuxControl_(other[jobs], copy));
  }
  Log::trace() << "ObsAuxControls<MODEL>::ObsAuxControls copy done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ObsAuxControls<MODEL>::~ObsAuxControls() {
  Log::trace() << "ObsAuxControls<MODEL>::~ObsAuxControls starting" << std::endl;
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) auxs_[jobs].reset();
  Log::trace() << "ObsAuxControls<MODEL>::~ObsAuxControls done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsAuxControls<MODEL>::read(const eckit::Configuration & conf) {
  Log::trace() << "ObsAuxControls<MODEL>::read starting" << std::endl;
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) auxs_[jobs]->read(conf);
  Log::trace() << "ObsAuxControls<MODEL>::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsAuxControls<MODEL>::write(const eckit::Configuration & conf) const {
  Log::trace() << "ObsAuxControls<MODEL>::write starting" << std::endl;
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) auxs_[jobs]->write(conf);
  Log::trace() << "ObsAuxControls<MODEL>::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double ObsAuxControls<MODEL>::norm() const {
  Log::trace() << "ObsAuxControls<MODEL>::norm starting" << std::endl;
  double zz =  static_cast<double>(0.0);
  std::size_t ii = 0;
  double norm;
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) {
    norm = auxs_[jobs]->norm();
    if (norm > 0.0) {
      zz += norm*norm;
      ++ii;
    }
  }
  Log::trace() << "ObsAuxControls<MODEL>::norm done" << std::endl;
  return std::sqrt(zz/ii);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsAuxControls<MODEL>::print(std::ostream & os) const {
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) os << *auxs_[jobs];
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSAUXCONTROLS_H_
