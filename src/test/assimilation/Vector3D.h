/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef TEST_ASSIMILATION_VECTOR3D_H_
#define TEST_ASSIMILATION_VECTOR3D_H_

namespace test {

  class Vector3D {
   public:
    Vector3D(const double& x,
             const double& y,
             const double& z);
    Vector3D& operator=(const Vector3D&);
    Vector3D& operator+=(const Vector3D&);
    Vector3D& operator-=(const Vector3D&);
    Vector3D& operator*=(const double);
    Vector3D& operator*=(const Vector3D&);
    Vector3D& operator/=(const Vector3D&);
    /// x -> x + mult * rhs
    void axpy(const double, const Vector3D&);
    double dot_product_with(const Vector3D&) const;
    double x() {return x_;}
    double y() {return y_;}
    double z() {return z_;}

   private:
    double x_;
    double y_;
    double z_;
  };
}  // namespace test

#endif  // TEST_ASSIMILATION_VECTOR3D_H_
