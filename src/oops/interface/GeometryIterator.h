/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_GEOMETRYITERATOR_H_
#define OOPS_INTERFACE_GEOMETRYITERATOR_H_

#include <iterator>
#include <memory>
#include <string>


#include "eckit/geometry/Point2.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

// -----------------------------------------------------------------------------
template<typename MODEL>
class GeometryIterator: public std::iterator<std::forward_iterator_tag,
                                             eckit::geometry::Point2>,
                        public util::Printable,
                        private util::ObjectCounter<GeometryIterator<MODEL>> {
  typedef typename MODEL::GeometryIterator GeometryIterator_;

 public:
  static const std::string classname() {return "oops::GeometryIterator";}

  GeometryIterator(const GeometryIterator&);
  explicit GeometryIterator(const GeometryIterator_&);
  ~GeometryIterator();

  bool operator==(const GeometryIterator&);
  bool operator!=(const GeometryIterator&);
  eckit::geometry::Point2 operator*();
  GeometryIterator operator++();

/// Interfacing
  const GeometryIterator_ & geometryiter() const {return *geometryiter_;}
  GeometryIterator_ & geometryiter() {return *geometryiter_; }

 private:
  void print(std::ostream &) const;
  std::unique_ptr<GeometryIterator_> geometryiter_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
GeometryIterator<MODEL>::GeometryIterator(const GeometryIterator& other) {
  Log::trace() << "GeometryIterator<MODEL>::GeometryIterator starting" << std::endl;
  util::Timer timer(classname(), "GeometryIterator");
  geometryiter_.reset(new GeometryIterator_(other.geometryiter()));
  Log::trace() << "GeometryIterator<MODEL>::GeometryIterator done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
GeometryIterator<MODEL>::GeometryIterator(const GeometryIterator_& iter) {
  Log::trace() << "GeometryIterator<MODEL>::GeometryIterator starting" << std::endl;
  util::Timer timer(classname(), "GeometryIterator");
  geometryiter_.reset(new GeometryIterator_(iter));
  Log::trace() << "GeometryIterator<MODEL>::GeometryIterator done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
GeometryIterator<MODEL>::~GeometryIterator() {
  Log::trace() << "GeometryIterator<MODEL>::~GeometryIterator starting" << std::endl;
  util::Timer timer(classname(), "~GeometryIterator");
  geometryiter_.reset();
  Log::trace() << "GeometryIterator<MODEL>::~GeometryIterator done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
bool GeometryIterator<MODEL>::operator==(const GeometryIterator& other) {
  Log::trace() << "GeometryIterator<MODEL>::operator== starting" << std::endl;
  util::Timer timer(classname(), "operator==");
  bool equals = (*geometryiter_ == other.geometryiter());
  Log::trace() << "GeometryIterator<MODEL>::operator== done" << std::endl;
  return equals;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
bool GeometryIterator<MODEL>::operator!=(const GeometryIterator& other) {
  Log::trace() << "GeometryIterator<MODEL>::operator!= starting" << std::endl;
  util::Timer timer(classname(), "operator!=");
  bool notequals = (*geometryiter_ != other.geometryiter());
  Log::trace() << "GeometryIterator<MODEL>::operator!= done" << std::endl;
  return notequals;
}


// -----------------------------------------------------------------------------

template<typename MODEL>
eckit::geometry::Point2 GeometryIterator<MODEL>::operator*() {
  Log::trace() << "GeometryIterator<MODEL>::operator* starting" << std::endl;
  util::Timer timer(classname(), "operator*");
  eckit::geometry::Point2 loc = *(*geometryiter_);
  Log::trace() << "GeometryIterator<MODEL>::operator* done" << std::endl;
  return loc;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
GeometryIterator<MODEL> GeometryIterator<MODEL>::operator++() {
  Log::trace() << "GeometryIterator<MODEL>::operator++ starting" << std::endl;
  util::Timer timer(classname(), "operator++");
  ++(*geometryiter_);
  Log::trace() << "GeometryIterator<MODEL>::operator++ done" << std::endl;
  return *this;
}


// -----------------------------------------------------------------------------

template<typename MODEL>
void GeometryIterator<MODEL>::print(std::ostream & os) const {
  Log::trace() << "GeometryIterator<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *geometryiter_;
  Log::trace() << "GeometryIterator<MODEL>::print done" << std::endl;
}


}  // namespace oops

#endif  // OOPS_INTERFACE_GEOMETRYITERATOR_H_
