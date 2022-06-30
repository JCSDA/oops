/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_UTIL_OBJECTCOUNTER_H_
#define OOPS_UTIL_OBJECTCOUNTER_H_

#include <memory>

#include "eckit/exception/Exceptions.h"

#include "oops/util/ObjectCountHelper.h"

namespace util {

// -----------------------------------------------------------------------------

template<typename T>
class ObjectCounter {
 public:
  ObjectCounter(): count_(ObjectCountHelper::create(T::classname())), bytes_(0)
    {count_->oneMore();}

  ObjectCounter(const ObjectCounter & other) : count_(other.count_), bytes_(0)
    {count_->oneMore();}

  ~ObjectCounter()
    {count_->oneLess(bytes_);}

 protected:
  size_t created() const {return count_->created();}

// Optionally set object size (in bytes)
  void setObjectSize(const size_t & bytes) {
    ASSERT(bytes_ == 0);  // can only be set once
    bytes_ = bytes;
    count_->setSize(bytes_);
  }

 private:
  std::shared_ptr<ObjectCountHelper> count_;
  size_t bytes_;
};

// -----------------------------------------------------------------------------

}  // namespace util

#endif  // OOPS_UTIL_OBJECTCOUNTER_H_
