/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_UTIL_DOT_PRODUCT_H_
#define OOPS_UTIL_DOT_PRODUCT_H_

/// Syntactic sugar to let us use a more mathematical notation for dot products.

template<class T>
inline
double dot_product(const T& x, const T& y) {
  return x.dot_product_with(y);
}

#endif  // OOPS_UTIL_DOT_PRODUCT_H_
