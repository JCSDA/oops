/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "oops/util/ObjectCountHelper.h"

#include <algorithm>
#include <iomanip>
#include <string>

#include "eckit/exception/Exceptions.h"
#include "oops/util/Logger.h"

namespace util {

// -----------------------------------------------------------------------------

std::map<std::string, std::shared_ptr<ObjectCountHelper> > ObjectCountHelper::counters_;

// -----------------------------------------------------------------------------

void ObjectCountHelper::start() {
  oops::Log::stats() << "ObjectCountHelper started." << std::endl;
}

// -----------------------------------------------------------------------------

void ObjectCountHelper::stop() {
  typedef std::map<std::string, std::shared_ptr<ObjectCountHelper> >::iterator it;
  oops::Log::stats() << " " << std::endl;
  oops::Log::stats() << "----------------------------------------------------------------------"
                     << "------------" << std::endl;
  oops::Log::stats() << "--------------------------- Object counts ----------------------------"
                     << "------------" << std::endl;
  oops::Log::stats() << "----------------------------------------------------------------------"
                     << "------------" << std::endl;
  oops::Log::stats() << std::setw(34) << std::left << " "
                     << std::setw(8) << std::right << "Total"
                     << std::setw(9) << std::right << "Simult."
                     << std::setw(7) << std::right << "Remain"
                     << std::setw(12) << std::right << "Avg (Mb)"
                     << std::setw(12) << std::right << "HWM (Mb)"
                     << std::endl;
  for (it jc = counters_.begin(); jc != counters_.end(); ++jc) {
    oops::Log::stats() << std::setw(32) << std::left << jc->first
                       << ": " << *(jc->second) << std::endl;
  }
  oops::Log::stats() << "----------------------------- Object counts --------------------------"
                     << "------------" << std::endl;
  counters_.clear();
}

// -----------------------------------------------------------------------------

std::shared_ptr<ObjectCountHelper> ObjectCountHelper::create(const std::string & cname) {
  std::shared_ptr<ObjectCountHelper> pcount;
  typedef std::map<std::string, std::shared_ptr<ObjectCountHelper> >::iterator it;
  it jj = counters_.find(cname);
  if (jj == counters_.end()) {
    pcount.reset(new ObjectCountHelper(cname));
    counters_[cname] = pcount;
  } else {
    pcount = jj->second;
  }
  return pcount;
}

// -----------------------------------------------------------------------------

ObjectCountHelper::ObjectCountHelper(const std::string & cname)
    : current_(0), created_(0), max_(0), bytes_(0), maxbytes_(0), totbytes_(0) {}

// -----------------------------------------------------------------------------

ObjectCountHelper::~ObjectCountHelper() {}

// -----------------------------------------------------------------------------

void ObjectCountHelper::oneMore() {
  ++current_;
  ++created_;
  max_ = std::max(max_, current_);
}

// -----------------------------------------------------------------------------

void ObjectCountHelper::oneLess(const size_t & bytes) {
  --current_;
  bytes_ -= bytes;
}

// -----------------------------------------------------------------------------

void ObjectCountHelper::setSize(const size_t & bytes) {
  bytes_ += bytes;
  maxbytes_ = std::max(maxbytes_, bytes_);
  totbytes_ += bytes;
}

// -----------------------------------------------------------------------------

void ObjectCountHelper::print(std::ostream & out) const {
  out << std::setw(8) << std::right << created_
      << std::setw(8) << std::right << max_;
  if (current_ > 0) {
    out << std::setw(8) << std::right << current_;
  } else {
    out << std::setw(8) << " ";
  }
  if (maxbytes_ > 0) {
    double size = static_cast<double>(totbytes_) / static_cast<double>(created_) / 1.e+6;
    out << std::setw(12) << std::right << std::fixed << std::setprecision(2) << size;
    double hwm = static_cast<double>(maxbytes_) / 1.e+6;
    out << std::setw(12) << std::right << std::fixed << std::setprecision(2) << hwm;
  }
}

// -----------------------------------------------------------------------------

}  // namespace util
