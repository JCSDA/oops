/*
 * (C) Copyright 2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_ASSIMILATION_VECTOR2D_H_
#define TEST_ASSIMILATION_VECTOR2D_H_

#include "oops/util/Printable.h"

namespace test {

  class Vector2D : public util::Printable {
   public:
    Vector2D(const double& x,
             const double& y);
    Vector2D(const Vector2D&, const bool copy = true);
    Vector2D& operator=(const Vector2D&);
    Vector2D& operator+=(const Vector2D&);
    Vector2D& operator-=(const Vector2D&);
    Vector2D& operator*=(const double);
    Vector2D& operator*=(const Vector2D&);
    Vector2D& operator/=(const Vector2D&);
    /// x -> x + mult * rhs
    void axpy(const double, const Vector2D&);
    double dot_product_with(const Vector2D&) const;
    void multiply(const Vector2D&, Vector2D&);
    double x() const {return x_;}
    double y() const {return y_;}

   private:
    void print(std::ostream & os) const;

    double x_;
    double y_;
  };
}  // namespace test

#endif  // TEST_ASSIMILATION_VECTOR2D_H_
