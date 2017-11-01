/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef UTIL_OBJECTCOUNTHELPER_H_
#define UTIL_OBJECTCOUNTHELPER_H_

#include <map>
#include <ostream>
#include <string>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include "util/Printable.h"

namespace util {

// -----------------------------------------------------------------------------

class ObjectCountHelper : public util::Printable,
                          private boost::noncopyable {
 public:
  static void start();
  static void stop();
  static boost::shared_ptr<ObjectCountHelper> create(const std::string &);

  ~ObjectCountHelper();
  void oneMore();
  void oneLess();

 private:
  static std::map< std::string, boost::shared_ptr<ObjectCountHelper> > * counters_;

  explicit ObjectCountHelper(const std::string &);
  void print(std::ostream &) const;

  unsigned int current_;
  unsigned int created_;
  unsigned int max_;
};

// -----------------------------------------------------------------------------

}  // namespace util

#endif  // UTIL_OBJECTCOUNTHELPER_H_
