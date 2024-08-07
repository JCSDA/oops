/*
 * (C) Copyright 2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_ASSIMILATION_VECTOR4D_H_
#define TEST_ASSIMILATION_VECTOR4D_H_

#include "oops/util/Printable.h"

namespace test {

  class Vector4D : public util::Printable {
   public:
    Vector4D(const double& x,
             const double& y,
             const double& z,
             const double& w);
    Vector4D(const Vector4D&, const bool copy = true);
    Vector4D& operator=(const Vector4D&);
    Vector4D& operator+=(const Vector4D&);
    Vector4D& operator-=(const Vector4D&);
    Vector4D& operator*=(const double);
    Vector4D& operator*=(const Vector4D&);
    Vector4D& operator/=(const Vector4D&);
    /// x -> x + mult * rhs
    void axpy(const double, const Vector4D&);
    double dot_product_with(const Vector4D&) const;
    void multiply(const Vector4D&, Vector4D&);
    double x() const {return x_;}
    double y() const {return y_;}
    double z() const {return z_;}
    double w() const {return w_;}

   private:
    void print(std::ostream & os) const;

    double x_;
    double y_;
    double z_;
    double w_;
  };
}  // namespace test

#endif  // TEST_ASSIMILATION_VECTOR4D_H_
