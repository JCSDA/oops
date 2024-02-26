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

#include "eckit/geometry/Point3.h"

#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

// -----------------------------------------------------------------------------
template<typename TRAIT>
class GeometryIterator: public util::Printable,
                        private util::ObjectCounter<GeometryIterator<TRAIT>> {
  typedef typename TRAIT::GeometryIterator GeometryIterator_;

 public:
  typedef eckit::geometry::Point3 value_type;
  typedef value_type& reference;
  typedef value_type* pointer;
  typedef std::forward_iterator_tag iterator_category;
  typedef std::ptrdiff_t difference_type;

  static const std::string classname() {return "oops::GeometryIterator";}

  GeometryIterator(const GeometryIterator&);
  explicit GeometryIterator(const GeometryIterator_&);


  ~GeometryIterator();

  bool operator==(const GeometryIterator&);
  bool operator!=(const GeometryIterator&);
  eckit::geometry::Point3 operator*() const;
  // pre-increment operator
  GeometryIterator& operator++();

/// Interfacing
  const GeometryIterator_ & geometryiter() const {return *geometryiter_;}

  GeometryIterator_ & geometryiter() {return *geometryiter_;}

 private:
  void print(std::ostream &) const;
  std::unique_ptr<GeometryIterator_> geometryiter_;
};

// -----------------------------------------------------------------------------

template<typename TRAIT>
GeometryIterator<TRAIT>::GeometryIterator(const GeometryIterator& other) {
  Log::trace() << "GeometryIterator<TRAIT>::GeometryIterator starting" << std::endl;
  util::Timer timer(classname(), "GeometryIterator");
  geometryiter_.reset(new GeometryIterator_(other.geometryiter()));
  Log::trace() << "GeometryIterator<TRAIT>::GeometryIterator done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename TRAIT>
GeometryIterator<TRAIT>::GeometryIterator(const GeometryIterator_& iter) {
  Log::trace() << "GeometryIterator<TRAIT>::GeometryIterator starting" << std::endl;
  util::Timer timer(classname(), "GeometryIterator");
  geometryiter_.reset(new GeometryIterator_(iter));
  Log::trace() << "GeometryIterator<TRAIT>::GeometryIterator done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename TRAIT>
GeometryIterator<TRAIT>::~GeometryIterator() {
  Log::trace() << "GeometryIterator<TRAIT>::~GeometryIterator starting" << std::endl;
  util::Timer timer(classname(), "~GeometryIterator");
  geometryiter_.reset();
  Log::trace() << "GeometryIterator<TRAIT>::~GeometryIterator done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename TRAIT>
bool GeometryIterator<TRAIT>::operator==(const GeometryIterator& other) {
  Log::trace() << "GeometryIterator<TRAIT>::operator== starting" << std::endl;
  util::Timer timer(classname(), "operator==");
  bool equals = (*geometryiter_ == other.geometryiter());
  Log::trace() << "GeometryIterator<TRAIT>::operator== done" << std::endl;
  return equals;
}

// -----------------------------------------------------------------------------

template<typename TRAIT>
bool GeometryIterator<TRAIT>::operator!=(const GeometryIterator& other) {
  Log::trace() << "GeometryIterator<TRAIT>::operator!= starting" << std::endl;
  util::Timer timer(classname(), "operator!=");
  bool notequals = (*geometryiter_ != other.geometryiter());
  Log::trace() << "GeometryIterator<TRAIT>::operator!= done" << std::endl;
  return notequals;
}


// -----------------------------------------------------------------------------

template<typename TRAIT>
eckit::geometry::Point3 GeometryIterator<TRAIT>::operator*() const {
  Log::trace() << "GeometryIterator<TRAIT>::operator* starting" << std::endl;
  util::Timer timer(classname(), "operator*");
  eckit::geometry::Point3 loc = *(*geometryiter_);
  Log::trace() << "GeometryIterator<TRAIT>::operator* done" << std::endl;
  return loc;
}

// -----------------------------------------------------------------------------

template<typename TRAIT>
GeometryIterator<TRAIT>& GeometryIterator<TRAIT>::operator++() {
  Log::trace() << "GeometryIterator<TRAIT>::operator++ starting" << std::endl;
  util::Timer timer(classname(), "operator++");
  ++(*geometryiter_);
  Log::trace() << "GeometryIterator<TRAIT>::operator++ done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename TRAIT>
void GeometryIterator<TRAIT>::print(std::ostream & os) const {
  Log::trace() << "GeometryIterator<TRAIT>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *geometryiter_;
  Log::trace() << "GeometryIterator<TRAIT>::print done" << std::endl;
}


}  // namespace oops

#endif  // OOPS_INTERFACE_GEOMETRYITERATOR_H_
