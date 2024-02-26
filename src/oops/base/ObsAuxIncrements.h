/*
 * (C) Copyright 2017-2019 UCAR
 * (C) Crown Copyright 2023, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSAUXINCREMENTS_H_
#define OOPS_BASE_OBSAUXINCREMENTS_H_

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsSpaces.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Holds a vector of ObsAuxIncrement
// -----------------------------------------------------------------------------

template <typename OBS>
class ObsAuxIncrements : public util::Printable,
                         public util::Serializable,
                         private util::ObjectCounter<ObsAuxIncrements<OBS> > {
  typedef ObsAuxIncrement<OBS>     ObsAuxIncrement_;
  typedef ObsAuxControls<OBS>      ObsAuxControls_;
  typedef ObsSpaces<OBS>           ObsSpaces_;

 public:
  static const std::string classname() {return "oops::ObsAuxIncrements";}

/// Constructor, destructor
  ObsAuxIncrements(const ObsSpaces_ &, const eckit::Configuration &);
  ObsAuxIncrements(const ObsAuxIncrements &, const bool copy = true);
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
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

 private:
  void print(std::ostream &) const override;
  std::vector<std::unique_ptr<ObsAuxIncrement_> > auxs_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsAuxControls<OBS> & operator+=(ObsAuxControls<OBS> & xx, const ObsAuxIncrements<OBS> & dx) {
  Log::trace() << "operator+=(ObsAuxControls, ObsAuxIncrements) starting" << std::endl;
  ASSERT(xx.size() == dx.size());
  for (std::size_t jobs = 0; jobs < xx.size(); ++jobs) {
    xx[jobs].obsauxcontrol() += dx[jobs].obsauxincrement();
  }
  Log::trace() << "operator+=(ObsAuxControls, ObsAuxIncrements) done" << std::endl;
  return xx;
}

// =============================================================================

template<typename OBS>
ObsAuxIncrements<OBS>::ObsAuxIncrements(const ObsSpaces_ & odb, const eckit::Configuration & conf)
  : auxs_(0)
{
  Log::trace() << "ObsAuxIncrements<OBS>::ObsAuxIncrements starting" << std::endl;
  size_t bytes = 0;
  std::vector<eckit::LocalConfiguration> obsconf = conf.getSubConfigurations();
  for (std::size_t jobs = 0; jobs < obsconf.size(); ++jobs) {
    eckit::LocalConfiguration obsauxconf = obsconf[jobs].getSubConfiguration("obs bias");
    auxs_.push_back(
      std::unique_ptr<ObsAuxIncrement_>(new ObsAuxIncrement_(odb[jobs], obsauxconf)));
    bytes += auxs_[jobs]->serialSize();
  }
  this->setObjectSize(bytes*sizeof(double));
  Log::trace() << "ObsAuxIncrements<OBS>::ObsAuxIncrements done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename OBS>
ObsAuxIncrements<OBS>::ObsAuxIncrements(const ObsAuxIncrements & other, const bool copy)
  : auxs_(other.size())
{
  Log::trace() << "ObsAuxIncrements<OBS>::ObsAuxIncrements copy starting" << std::endl;
  size_t bytes = 0;
  ASSERT(size() == other.size());
  for (std::size_t jobs = 0; jobs < other.size(); ++jobs) {
    auxs_[jobs].reset(new ObsAuxIncrement_(other[jobs], copy));
    bytes += auxs_[jobs]->serialSize();
  }
  this->setObjectSize(bytes*sizeof(double));
  Log::trace() << "ObsAuxIncrements<OBS>::ObsAuxIncrements copy done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename OBS>
ObsAuxIncrements<OBS>::~ObsAuxIncrements() {
  Log::trace() << "ObsAuxIncrements<OBS>::~ObsAuxIncrements starting" << std::endl;
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) auxs_[jobs].reset();
  Log::trace() << "ObsAuxIncrements<OBS>::~ObsAuxIncrements done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename OBS>
void ObsAuxIncrements<OBS>::diff(const ObsAuxControls_ & x1, const ObsAuxControls_ & x2) {
  Log::trace() << "ObsAuxIncrements<OBS>::diff starting" << std::endl;
  ASSERT(x1.size() == x2.size() &&  size() == x2.size());
  for (std::size_t jobs = 0; jobs < size(); ++jobs) {
    auxs_[jobs]->diff(x1[jobs], x2[jobs]);
  }
  Log::trace() << "ObsAuxIncrements<OBS>::diff done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename OBS>
void ObsAuxIncrements<OBS>::zero() {
  Log::trace() << "ObsAuxIncrements<OBS>::zero starting" << std::endl;
  for (std::size_t jobs = 0; jobs < size(); ++jobs) {
    auxs_[jobs]->zero();
  }
  Log::trace() << "ObsAuxIncrements<OBS>::zero done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename OBS>
ObsAuxIncrements<OBS> & ObsAuxIncrements<OBS>::operator=(const ObsAuxIncrements & rhs) {
  Log::trace() << "ObsAuxIncrements<OBS>::operator= starting" << std::endl;
  ASSERT(size() == rhs.size());
  for (std::size_t jobs = 0; jobs < size(); ++jobs) {
    *auxs_[jobs] = rhs[jobs];
  }
  Log::trace() << "ObsAuxIncrements<OBS>::operator= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename OBS>
ObsAuxIncrements<OBS> & ObsAuxIncrements<OBS>::operator+=(const ObsAuxIncrements & rhs) {
  Log::trace() << "ObsAuxIncrements<OBS>::operator+= starting" << std::endl;
  ASSERT(size() == rhs.size());
  for (std::size_t jobs = 0; jobs < size(); ++jobs) {
    *auxs_[jobs] += rhs[jobs];
  }
  Log::trace() << "ObsAuxIncrements<OBS>::operator+= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename OBS>
ObsAuxIncrements<OBS> & ObsAuxIncrements<OBS>::operator-=(const ObsAuxIncrements & rhs) {
  Log::trace() << "ObsAuxIncrements<OBS>::operator-= starting" << std::endl;
  ASSERT(size() == rhs.size());
  for (std::size_t jobs = 0; jobs < size(); ++jobs) {
    *auxs_[jobs] -= rhs[jobs];
  }
  Log::trace() << "ObsAuxIncrements<OBS>::operator-= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename OBS>
ObsAuxIncrements<OBS> & ObsAuxIncrements<OBS>::operator*=(const double & zz) {
  Log::trace() << "ObsAuxIncrements<OBS>::operator*= starting" << std::endl;
  for (std::size_t jobs = 0; jobs < size(); ++jobs) {
    *auxs_[jobs] *= zz;
  }
  Log::trace() << "ObsAuxIncrements<OBS>::operator*= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename OBS>
void ObsAuxIncrements<OBS>::axpy(const double & zz, const ObsAuxIncrements & dx) {
  Log::trace() << "ObsAuxIncrements<OBS>::axpy starting" << std::endl;
  ASSERT(size() == dx.size());
  for (std::size_t jobs = 0; jobs < size(); ++jobs) {
    auxs_[jobs]->axpy(zz, dx[jobs]);
  }
  Log::trace() << "ObsAuxIncrements<OBS>::axpy done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename OBS>
double ObsAuxIncrements<OBS>::dot_product_with(const ObsAuxIncrements & dx) const {
  Log::trace() << "ObsAuxIncrements<OBS>::dot_product_with starting" << std::endl;
  ASSERT(size() == dx.size());
  double zz = static_cast<double>(0.0);
  for (std::size_t jobs = 0; jobs < size(); ++jobs) {
    zz += auxs_[jobs]->dot_product_with(dx[jobs]);
  }
  Log::trace() << "ObsAuxIncrements<OBS>::dot_product_with done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
template<typename OBS>
void ObsAuxIncrements<OBS>::read(const eckit::Configuration & conf) {
  Log::trace() << "ObsAuxIncrements<OBS>::read starting" << std::endl;
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) auxs_[jobs]->read(conf);
  Log::trace() << "ObsAuxIncrements<OBS>::read done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename OBS>
void ObsAuxIncrements<OBS>::write(const eckit::Configuration & conf) const {
  Log::trace() << "ObsAuxIncrements<OBS>::write starting" << std::endl;
  const bool write = conf.getBool("write increment", false);
  if (write) {
    for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) {
      auxs_[jobs]->write(conf);
    }
  }
  Log::trace() << "ObsAuxIncrements<OBS>::write done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename OBS>
double ObsAuxIncrements<OBS>::norm() const {
  Log::trace() << "ObsAuxIncrements<OBS>::norm starting" << std::endl;
  double zz = static_cast<double>(0.0);
  for (std::size_t jobs = 0; jobs < size(); ++jobs) {
    zz += auxs_[jobs]->norm();
  }
  Log::trace() << "ObsAuxIncrements<OBS>::norm done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
template<typename OBS>
size_t ObsAuxIncrements<OBS>::serialSize() const {
  Log::trace() << "ObsAuxIncrements<OBS>::serialSize starting" << std::endl;
  size_t ss = 0;
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) {
    ss += auxs_[jobs]->serialSize();
  }
  Log::trace() << "ObsAuxIncrements<OBS>::serialSize done" << std::endl;
  return ss;
}
// -----------------------------------------------------------------------------
template<typename OBS>
void ObsAuxIncrements<OBS>::serialize(std::vector<double> & vect) const {
  Log::trace() << "ObsAuxIncrements<OBS>::serialize starting" << std::endl;
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) auxs_[jobs]->serialize(vect);
  Log::trace() << "ObsAuxIncrements<OBS>::serialize done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename OBS>
void ObsAuxIncrements<OBS>::deserialize(const std::vector<double> & vect, size_t & index) {
  Log::trace() << "ObsAuxIncrements<OBS>::deserialize starting" << std::endl;
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) {
    auxs_[jobs]->deserialize(vect, index);
  }
  Log::trace() << "ObsAuxIncrements<OBS>::deserialize done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename OBS>
void ObsAuxIncrements<OBS>::print(std::ostream & os) const {
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) os << *auxs_[jobs];
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSAUXINCREMENTS_H_
