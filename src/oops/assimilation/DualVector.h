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

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/base/Departures.h"
#include "oops/base/GeneralizedDepartures.h"
#include "oops/interface/Increment.h"
#include "util/dot_product.h"

namespace oops {

/// Container of dual space vectors for all terms of the cost function.
/*!
 * Contains dual space vectors for all terms of the cost function,
 * that is Departures for Jo, an Increment_ for Jc, a ControlIncrement
 * for Jb and Jq.
 */

// -----------------------------------------------------------------------------
template<typename MODEL> class DualVector {
  typedef Increment<MODEL>           Increment_;
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef Departures<MODEL>          Departures_;

 public:
  DualVector(): dxjb_(), dxjo_(), dxjc_(), ijo_(), ijc_(), size_(0) {}
  DualVector(const DualVector &);
  ~DualVector() {}

// Access increment
  void dx(CtrlInc_ * dx) {dxjb_.reset(dx);}
  const CtrlInc_ & dx() const {return *dxjb_;}
  CtrlInc_ & dx() {return *dxjb_;}

// Store and retrieve other elements (takes ownership)
  void append(GeneralizedDepartures *);
  boost::shared_ptr<const GeneralizedDepartures> getv(const unsigned) const;

// Linear algebra
  DualVector & operator=(const DualVector &);
  DualVector & operator+=(const DualVector &);
  DualVector & operator-=(const DualVector &);
  DualVector & operator*=(const double);
  void zero();
  void axpy(const double, const DualVector &);
  double dot_product_with(const DualVector &) const;

// Clear everything
  void clear();

 private:
  bool compatible(const DualVector & other) const;

  boost::scoped_ptr<CtrlInc_>   dxjb_;
  std::vector<boost::shared_ptr<Departures_> >    dxjo_;
  std::vector<boost::shared_ptr<Increment_> > dxjc_;
  std::vector<unsigned> ijo_;
  std::vector<unsigned> ijc_;
  unsigned size_;
};

// =============================================================================

template<typename MODEL>
DualVector<MODEL>::DualVector(const DualVector & other)
  : dxjb_(), dxjo_(), dxjc_(),
    ijo_(other.ijo_), ijc_(other.ijc_), size_(other.size_)
{
  if (other.dxjb_ != 0) {
    dxjb_.reset(new CtrlInc_(*other.dxjb_));
  }
  for (unsigned jj = 0; jj < other.dxjo_.size(); ++jj) {
    boost::shared_ptr<Departures_> pd(new Departures_(*other.dxjo_[jj]));
    dxjo_.push_back(pd);
  }
  for (unsigned jj = 0; jj < other.dxjc_.size(); ++jj) {
    boost::shared_ptr<Increment_> pi(new Increment_(*other.dxjc_[jj]));
    dxjc_.push_back(pi);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void DualVector<MODEL>::clear() {
  dxjb_.reset();
  dxjo_.clear();
  dxjc_.clear();
  ijo_.clear();
  ijc_.clear();
  size_ = 0;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void DualVector<MODEL>::append(GeneralizedDepartures * pv) {
// Since there is no duck-typing in C++, we do it manually.
  Increment_ * pi = dynamic_cast<Increment_*>(pv);
  if (pi != 0) {
    boost::shared_ptr<Increment_> si(pi);
    dxjc_.push_back(si);
    ijc_.push_back(size_);
  }
  Departures_ * pd = dynamic_cast<Departures_*>(pv);
  if (pd != 0) {
    boost::shared_ptr<Departures_> sd(pd);
    dxjo_.push_back(sd);
    ijo_.push_back(size_);
  }
  ASSERT(pi != 0 || pd != 0);
  ++size_;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
boost::shared_ptr<const GeneralizedDepartures>
DualVector<MODEL>::getv(const unsigned ii) const {
  ASSERT(ii < size_);
  boost::shared_ptr<const GeneralizedDepartures> pv;
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
template<typename MODEL>
DualVector<MODEL> & DualVector<MODEL>::operator=(const DualVector & rhs) {
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
template<typename MODEL>
DualVector<MODEL> & DualVector<MODEL>::operator+=(const DualVector & rhs) {
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
template<typename MODEL>
DualVector<MODEL> & DualVector<MODEL>::operator-=(const DualVector & rhs) {
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
template<typename MODEL>
DualVector<MODEL> & DualVector<MODEL>::operator*=(const double zz) {
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
template<typename MODEL>
void DualVector<MODEL>::zero() {
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
template<typename MODEL>
void DualVector<MODEL>::axpy(const double zz, const DualVector & rhs) {
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
template<typename MODEL>
double DualVector<MODEL>::dot_product_with(const DualVector & x2) const {
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
template<typename MODEL>
bool DualVector<MODEL>::compatible(const DualVector & other) const {
  bool lcheck = (dxjb_ == 0) == (other.dxjb_ == 0)
             && (dxjo_.size() == other.dxjo_.size())
             && (dxjc_.size() == other.dxjc_.size());
  return lcheck;
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_ASSIMILATION_DUALVECTOR_H_
