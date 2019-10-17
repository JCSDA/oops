/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSAUXINCREMENTS_H_
#define OOPS_BASE_OBSAUXINCREMENTS_H_

#include <iostream>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "oops/base/ObsAuxControls.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsAuxIncrements : public util::Printable {
  typedef ObsAuxIncrement<MODEL>               ObsAuxIncrement_;
  typedef ObsAuxControls<MODEL>                ObsAuxControls_;

 public:
  static const std::string classname() {return "oops::ObsAuxIncrements";}

/// Constructor, destructor
  explicit ObsAuxIncrements(const eckit::Configuration &);
  ObsAuxIncrements(const ObsAuxIncrements &, const bool copy = true);
  ObsAuxIncrements(const ObsAuxIncrements &, const eckit::Configuration &);
  ~ObsAuxIncrements();

/// Access
  std::size_t size() const {return auxs_.size();}
  const ObsAuxIncrement_ & operator[](const std::size_t ii) const {return *auxs_.at(ii);}
  ObsAuxIncrement_ & operator[](const std::size_t ii) {return *auxs_.at(ii);}

/// Linear algebra operators
  void diff(const ObsAuxControls_ &, const ObsAuxControls_ &);
  void zero();
  ObsAuxIncrements & operator=(const ObsAuxIncrements &);
  ObsAuxIncrements & operator+=(const ObsAuxIncrements &);
  ObsAuxIncrements & operator-=(const ObsAuxIncrements &);
  ObsAuxIncrements & operator*=(const double &);
  void axpy(const double &, const ObsAuxIncrements &);
  double dot_product_with(const ObsAuxIncrements &) const;

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;

/// Serialize-Deserialize
  size_t serialSize() const;
  void serialize(std::vector<double> &) const;
  void deserialize(const std::vector<double>, size_t &);

 private:
  void print(std::ostream &) const;
  std::vector<boost::shared_ptr<ObsAuxIncrement_> > auxs_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsAuxControls<MODEL> & operator+=(ObsAuxControls<MODEL> & xx, const ObsAuxIncrements<MODEL> & dx) {
  Log::trace() << "operator+=(ObsAuxControls, ObsAuxIncrements) starting" << std::endl;
  ASSERT(xx.size() == dx.size());
  for (std::size_t jobs = 0; jobs < xx.size(); ++jobs) {
    xx[jobs].obsauxcontrol() += dx[jobs].obsauxincrement();
  }
  Log::trace() << "operator+=(ObsAuxControls, ObsAuxIncrements) done" << std::endl;
  return xx;
}

// =============================================================================

template<typename MODEL>
ObsAuxIncrements<MODEL>::ObsAuxIncrements(const eckit::Configuration & conf)
  : auxs_(0)
{
  std::vector<eckit::LocalConfiguration> obsconf;
  conf.get("ObsTypes", obsconf);
  for (std::size_t jobs = 0; jobs < obsconf.size(); ++jobs) {
    boost::shared_ptr<ObsAuxIncrement_>
          tmp(new ObsAuxIncrement_(obsconf[jobs]));
    auxs_.push_back(tmp);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ObsAuxIncrements<MODEL>::ObsAuxIncrements(const ObsAuxIncrements & other, const bool copy)
  : auxs_(other.size())
{
  Log::trace() << "ObsAuxIncrements<MODEL>::ObsAuxIncrements copy starting" << std::endl;
  ASSERT(size() == other.size());
  for (std::size_t jobs = 0; jobs < other.size(); ++jobs) {
    auxs_[jobs].reset(new ObsAuxIncrement_(other[jobs], copy));
  }
  Log::trace() << "ObsAuxIncrements<MODEL>::ObsAuxIncrements copy done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ObsAuxIncrements<MODEL>::ObsAuxIncrements(const ObsAuxIncrements & other,
                                        const eckit::Configuration & conf) : auxs_(other.size())
{
  Log::trace() << "ObsAuxIncrements<MODEL>::ObsAuxIncrements interpolated starting" << std::endl;
  std::vector<eckit::LocalConfiguration> obsconf;
  ASSERT(size() == other.size());
  for (std::size_t jobs = 0; jobs < other.size(); ++jobs) {
    auxs_[jobs].reset(new ObsAuxIncrement_(other[jobs]));
  }
  Log::trace() << "ObsAuxIncrements<MODEL>::ObsAuxIncrements interpolated done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ObsAuxIncrements<MODEL>::~ObsAuxIncrements() {
  Log::trace() << "ObsAuxIncrements<MODEL>::~ObsAuxIncrements starting" << std::endl;
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) auxs_[jobs].reset();
  Log::trace() << "ObsAuxIncrements<MODEL>::~ObsAuxIncrements done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ObsAuxIncrements<MODEL>::diff(const ObsAuxControls_ & x1, const ObsAuxControls_ & x2) {
  Log::trace() << "ObsAuxIncrements<MODEL>::diff starting" << std::endl;
  ASSERT(x1.size() == x2.size() &&  size() == x2.size());
  for (std::size_t jobs = 0; jobs < size(); ++jobs) {
    auxs_[jobs]->diff(x1[jobs], x2[jobs]);
  }
  Log::trace() << "ObsAuxIncrements<MODEL>::diff done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ObsAuxIncrements<MODEL>::zero() {
  Log::trace() << "ObsAuxIncrements<MODEL>::zero starting" << std::endl;
  for (std::size_t jobs = 0; jobs < size(); ++jobs) {
    auxs_[jobs]->zero();
  }
  Log::trace() << "ObsAuxIncrements<MODEL>::zero done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ObsAuxIncrements<MODEL> & ObsAuxIncrements<MODEL>::operator=(const ObsAuxIncrements & rhs) {
  Log::trace() << "ObsAuxIncrements<MODEL>::operator= starting" << std::endl;
  ASSERT(size() == rhs.size());
  for (std::size_t jobs = 0; jobs < size(); ++jobs) {
    *auxs_[jobs] = rhs[jobs];
  }
  Log::trace() << "ObsAuxIncrements<MODEL>::operator= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ObsAuxIncrements<MODEL> & ObsAuxIncrements<MODEL>::operator+=(const ObsAuxIncrements & rhs) {
  Log::trace() << "ObsAuxIncrements<MODEL>::operator+= starting" << std::endl;
  ASSERT(size() == rhs.size());
  for (std::size_t jobs = 0; jobs < size(); ++jobs) {
    *auxs_[jobs] += rhs[jobs];
  }
  Log::trace() << "ObsAuxIncrements<MODEL>::operator+= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ObsAuxIncrements<MODEL> & ObsAuxIncrements<MODEL>::operator-=(const ObsAuxIncrements & rhs) {
  Log::trace() << "ObsAuxIncrements<MODEL>::operator-= starting" << std::endl;
  ASSERT(size() == rhs.size());
  for (std::size_t jobs = 0; jobs < size(); ++jobs) {
    *auxs_[jobs] -= rhs[jobs];
  }
  Log::trace() << "ObsAuxIncrements<MODEL>::operator-= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ObsAuxIncrements<MODEL> & ObsAuxIncrements<MODEL>::operator*=(const double & zz) {
  Log::trace() << "ObsAuxIncrements<MODEL>::operator*= starting" << std::endl;
  for (std::size_t jobs = 0; jobs < size(); ++jobs) {
    *auxs_[jobs] *= zz;
  }
  Log::trace() << "ObsAuxIncrements<MODEL>::operator*= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ObsAuxIncrements<MODEL>::axpy(const double & zz, const ObsAuxIncrements & dx) {
  Log::trace() << "ObsAuxIncrements<MODEL>::axpy starting" << std::endl;
  ASSERT(size() == dx.size());
  for (std::size_t jobs = 0; jobs < size(); ++jobs) {
    auxs_[jobs]->axpy(zz, dx[jobs]);
  }
  Log::trace() << "ObsAuxIncrements<MODEL>::axpy done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
double ObsAuxIncrements<MODEL>::dot_product_with(const ObsAuxIncrements & dx) const {
  Log::trace() << "ObsAuxIncrements<MODEL>::dot_product_with starting" << std::endl;
  ASSERT(size() == dx.size());
  double zz = static_cast<double>(0.0);
  for (std::size_t jobs = 0; jobs < size(); ++jobs) {
    zz += auxs_[jobs]->dot_product_with(dx[jobs]);
  }
  Log::trace() << "ObsAuxIncrements<MODEL>::dot_product_with done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ObsAuxIncrements<MODEL>::read(const eckit::Configuration & conf) {
  Log::trace() << "ObsAuxIncrements<MODEL>::read starting" << std::endl;
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) auxs_[jobs]->read(conf);
  Log::trace() << "ObsAuxIncrements<MODEL>::read done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ObsAuxIncrements<MODEL>::write(const eckit::Configuration & conf) const {
  Log::trace() << "ObsAuxIncrements<MODEL>::write starting" << std::endl;
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) auxs_[jobs]->write(conf);
  Log::trace() << "ObsAuxIncrements<MODEL>::write done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
double ObsAuxIncrements<MODEL>::norm() const {
  Log::trace() << "ObsAuxIncrements<MODEL>::norm starting" << std::endl;
  double zz = static_cast<double>(0.0);
  for (std::size_t jobs = 0; jobs < size(); ++jobs) {
    zz += auxs_[jobs]->norm();
  }
  Log::trace() << "ObsAuxIncrements<MODEL>::norm done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
size_t ObsAuxIncrements<MODEL>::serialSize() const {
  Log::trace() << "ObsAuxIncrements<MODEL>::serialSize starting" << std::endl;
  size_t ss = 0;
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) {
    ss += auxs_[jobs]->serialSize();
  }
  Log::trace() << "ObsAuxIncrements<MODEL>::serialSize done" << std::endl;
  return ss;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ObsAuxIncrements<MODEL>::serialize(std::vector<double> & vect) const {
  Log::trace() << "ObsAuxIncrements<MODEL>::serialize starting" << std::endl;
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) auxs_[jobs]->serialize(vect);
  Log::trace() << "ObsAuxIncrements<MODEL>::serialize done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ObsAuxIncrements<MODEL>::deserialize(const std::vector<double> vect, size_t & index) {
  Log::trace() << "ObsAuxIncrements<MODEL>::deserialize starting" << std::endl;
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) {
    auxs_[jobs]->deserialize(vect, index);
  }
  Log::trace() << "ObsAuxIncrements<MODEL>::deserialize done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL>
void ObsAuxIncrements<MODEL>::print(std::ostream & os) const {
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) os << *auxs_[jobs];
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSAUXINCREMENTS_H_
