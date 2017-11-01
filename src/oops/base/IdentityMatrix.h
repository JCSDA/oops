/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_IDENTITYMATRIX_H_
#define OOPS_BASE_IDENTITYMATRIX_H_

namespace oops {

/// Identity matrix
template<typename VECTOR>
class IdentityMatrix {
 public:
  void multiply(const VECTOR & a, VECTOR & b) const {b = a;}
};

}  // namespace oops

#endif  // OOPS_BASE_IDENTITYMATRIX_H_
