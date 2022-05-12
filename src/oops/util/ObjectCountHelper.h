/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_UTIL_OBJECTCOUNTHELPER_H_
#define OOPS_UTIL_OBJECTCOUNTHELPER_H_

#include <map>
#include <memory>
#include <ostream>
#include <string>

#include <boost/noncopyable.hpp>
#include "oops/util/Printable.h"

namespace util {

// -----------------------------------------------------------------------------

class ObjectCountHelper : public util::Printable,
                          private boost::noncopyable {
 public:
  static void start();
  static void stop();
  static std::shared_ptr<ObjectCountHelper> create(const std::string &);

  ~ObjectCountHelper();
  void oneMore();
  void oneLess(const size_t &);
  void setSize(const size_t &);
  size_t created() const {return created_;}

 private:
  static std::map< std::string, std::shared_ptr<ObjectCountHelper> > counters_;

  explicit ObjectCountHelper(const std::string &);
  void print(std::ostream &) const;

  size_t current_;
  size_t created_;
  size_t max_;
  size_t bytes_;
  size_t maxbytes_;
  size_t totbytes_;
};

// -----------------------------------------------------------------------------

}  // namespace util

#endif  // OOPS_UTIL_OBJECTCOUNTHELPER_H_
