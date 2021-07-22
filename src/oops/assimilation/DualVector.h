/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_DUALVECTOR_H_
#define OOPS_ASSIMILATION_DUALVECTOR_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/base/Departures.h"
#include "oops/base/GeneralizedDepartures.h"
#include "oops/base/Increment.h"
#include "oops/util/dot_product.h"

namespace oops {

/// Container of dual space vectors for all terms of the cost function.
/*!
 * Contains dual space vectors for all terms of the cost function,
 * that is Departures for Jo, an Increment_ for Jc, a ControlIncrement
 * for Jb and Jq.
 */

// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS> class DualVector {
  typedef Increment<MODEL>                Increment_;
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef Departures<OBS>                 Departures_;

 public:
  DualVector(): dxjb_(), dxjo_(), dxjc_(), ijo_(), ijc_(), size_(0) {}
  DualVector(const DualVector &);
  ~DualVector() {}

// Access increment
  void dx(CtrlInc_ * dx) {dxjb_.reset(dx);}
  const CtrlInc_ & dx() const {return *dxjb_;}
  CtrlInc_ & dx() {return *dxjb_;}

// Store and retrieve other elements (takes ownership)
  void append(std::unique_ptr<GeneralizedDepartures> &&);
  std::shared_ptr<const GeneralizedDepartures> getv(const unsigned) const;

// Linear algebra
  DualVector & operator=(const DualVector &);
  DualVector & operator+=(const DualVector &);
  DualVector & operator-=(const DualVector &);
  DualVector & operator*=(const double);
  void zero();
  void axpy(const double, const DualVector &);
  double dot_product_with(const DualVector &) const;

  void saveDep(const std::string &) const;

// Clear everything
  void clear();

 private:
  bool compatible(const DualVector & other) const;

  std::unique_ptr<CtrlInc_>   dxjb_;
  std::vector<std::shared_ptr<Departures_> >    dxjo_;
  std::vector<std::shared_ptr<Increment_> > dxjc_;
  std::vector<unsigned> ijo_;
  std::vector<unsigned> ijc_;
  unsigned size_;
};

// =============================================================================

template<typename MODEL, typename OBS>
DualVector<MODEL, OBS>::DualVector(const DualVector & other)
  : dxjb_(), dxjo_(), dxjc_(),
    ijo_(other.ijo_), ijc_(other.ijc_), size_(other.size_)
{
  if (other.dxjb_ != 0) {
    dxjb_.reset(new CtrlInc_(*other.dxjb_));
  }
  for (unsigned jj = 0; jj < other.dxjo_.size(); ++jj) {
    std::shared_ptr<Departures_> pd(new Departures_(*other.dxjo_[jj]));
    dxjo_.push_back(pd);
  }
  for (unsigned jj = 0; jj < other.dxjc_.size(); ++jj) {
    std::shared_ptr<Increment_> pi(new Increment_(*other.dxjc_[jj]));
    dxjc_.push_back(pi);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
void DualVector<MODEL, OBS>::clear() {
  dxjb_.reset();
  dxjo_.clear();
  dxjc_.clear();
  ijo_.clear();
  ijc_.clear();
  size_ = 0;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
void DualVector<MODEL, OBS>::append(std::unique_ptr<GeneralizedDepartures> && uv) {
// Since there is no duck-typing in C++, we do it manually.
  std::shared_ptr<GeneralizedDepartures> sv = std::move(uv);
  std::shared_ptr<Increment_> si = std::dynamic_pointer_cast<Increment_>(sv);
  if (si != nullptr) {
    dxjc_.push_back(si);
    ijc_.push_back(size_);
  }
  std::shared_ptr<Departures_> sd = std::dynamic_pointer_cast<Departures_>(sv);
  if (sd != nullptr) {
    dxjo_.push_back(sd);
    ijo_.push_back(size_);
  }
  ASSERT(si != nullptr || sd != nullptr);
  ++size_;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
std::shared_ptr<const GeneralizedDepartures>
DualVector<MODEL, OBS>::getv(const unsigned ii) const {
  ASSERT(ii < size_);
  std::shared_ptr<const GeneralizedDepartures> pv;
  for (unsigned jj = 0; jj < ijo_.size(); ++jj) {
    if (ijo_[jj] == ii) pv = dxjo_[jj];
  }
  for (unsigned jj = 0; jj < ijc_.size(); ++jj) {
    if (ijc_[jj] == ii) pv = dxjc_[jj];
  }
  ASSERT(pv != 0);
  return pv;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
DualVector<MODEL, OBS> & DualVector<MODEL, OBS>::operator=(const DualVector & rhs) {
  ASSERT(this->compatible(rhs));
  if (dxjb_ != 0) {
    *dxjb_ = *rhs.dxjb_;
  }
  for (unsigned jj = 0; jj < dxjo_.size(); ++jj) {
    *dxjo_[jj] = *rhs.dxjo_[jj];
  }
  for (unsigned jj = 0; jj < dxjc_.size(); ++jj) {
    *dxjc_[jj] = *rhs.dxjc_[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
DualVector<MODEL, OBS> & DualVector<MODEL, OBS>::operator+=(const DualVector & rhs) {
  ASSERT(this->compatible(rhs));
  if (dxjb_ != 0) {
    *dxjb_ += *rhs.dxjb_;
  }
  for (unsigned jj = 0; jj < dxjo_.size(); ++jj) {
    *dxjo_[jj] += *rhs.dxjo_[jj];
  }
  for (unsigned jj = 0; jj < dxjc_.size(); ++jj) {
    *dxjc_[jj] += *rhs.dxjc_[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
DualVector<MODEL, OBS> & DualVector<MODEL, OBS>::operator-=(const DualVector & rhs) {
  ASSERT(this->compatible(rhs));
  if (dxjb_ != 0) {
    *dxjb_ -= *rhs.dxjb_;
  }
  for (unsigned jj = 0; jj < dxjo_.size(); ++jj) {
    *dxjo_[jj] -= *rhs.dxjo_[jj];
  }
  for (unsigned jj = 0; jj < dxjc_.size(); ++jj) {
    *dxjc_[jj] -= *rhs.dxjc_[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
DualVector<MODEL, OBS> & DualVector<MODEL, OBS>::operator*=(const double zz) {
  if (dxjb_ != 0) {
    *dxjb_ *= zz;
  }
  for (unsigned jj = 0; jj < dxjo_.size(); ++jj) {
    *dxjo_[jj] *= zz;
  }
  for (unsigned jj = 0; jj < dxjc_.size(); ++jj) {
    *dxjc_[jj] *= zz;
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
void DualVector<MODEL, OBS>::zero() {
  if (dxjb_ != 0) {
    dxjb_->zero();
  }
  for (unsigned jj = 0; jj < dxjo_.size(); ++jj) {
    dxjo_[jj]->zero();
  }
  for (unsigned jj = 0; jj < dxjc_.size(); ++jj) {
    dxjc_[jj]->zero();
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
void DualVector<MODEL, OBS>::axpy(const double zz, const DualVector & rhs) {
  ASSERT(this->compatible(rhs));
  if (dxjb_ != 0) {
    dxjb_->axpy(zz, *rhs.dxjb_);
  }
  for (unsigned jj = 0; jj < dxjo_.size(); ++jj) {
    dxjo_[jj]->axpy(zz, *rhs.dxjo_[jj]);
  }
  for (unsigned jj = 0; jj < dxjc_.size(); ++jj) {
    dxjc_[jj]->axpy(zz, *rhs.dxjc_[jj]);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
double DualVector<MODEL, OBS>::dot_product_with(const DualVector & x2) const {
  ASSERT(this->compatible(x2));
  double zz = 0.0;
  if (dxjb_ != 0) {
    zz += dot_product(*dxjb_, *x2.dxjb_);
  }
  for (unsigned jj = 0; jj < dxjo_.size(); ++jj) {
    zz += dot_product(*dxjo_[jj], *x2.dxjo_[jj]);
  }
  for (unsigned jj = 0; jj < dxjc_.size(); ++jj) {
    zz += dot_product(*dxjc_[jj], *x2.dxjc_[jj]);
  }
  return zz;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
bool DualVector<MODEL, OBS>::compatible(const DualVector & other) const {
  bool lcheck = (dxjb_ == 0) == (other.dxjb_ == 0)
             && (dxjo_.size() == other.dxjo_.size())
             && (dxjc_.size() == other.dxjc_.size());
  return lcheck;
}
// -----------------------------------------------------------------------------
template<typename MODEL, typename OBS>
void DualVector<MODEL, OBS>::saveDep(const std::string & name) const {
  for (unsigned jj = 0; jj < dxjo_.size(); ++jj) {
    dxjo_[jj]->save(name);
  }
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_ASSIMILATION_DUALVECTOR_H_
