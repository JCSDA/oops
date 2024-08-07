/*
 * (C) Copyright 2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "test/assimilation/Vector4D.h"

namespace test {

  Vector4D::Vector4D(const double& x,
                     const double& y,
                     const double& z,
                     const double& w)
    : x_(x), y_(y), z_(z), w_(w)
  {}

  Vector4D::Vector4D(const Vector4D& rhs,
                     const bool copy)
    : x_(rhs.x_), y_(rhs.y_), z_(rhs.z_), w_(rhs.w_)
  {
    // NB the argument 'copy' is present in order to match the signature
    // required by various minimiser functions.
  }

  Vector4D& Vector4D::operator=(const Vector4D& rhs)
  {
    this->x_ = rhs.x_;
    this->y_ = rhs.y_;
    this->z_ = rhs.z_;
    this->w_ = rhs.w_;
    return *this;
  }

  Vector4D& Vector4D::operator+=(const Vector4D& rhs)
  {
    this->x_ += rhs.x_;
    this->y_ += rhs.y_;
    this->z_ += rhs.z_;
    this->w_ += rhs.w_;
    return *this;
  }

  Vector4D& Vector4D::operator-=(const Vector4D& rhs)
  {
    this->x_ -= rhs.x_;
    this->y_ -= rhs.y_;
    this->z_ -= rhs.z_;
    this->w_ -= rhs.w_;
    return *this;
  }

  Vector4D& Vector4D::operator*=(const double mult)
  {
    this->x_ *= mult;
    this->y_ *= mult;
    this->z_ *= mult;
    this->w_ *= mult;
    return *this;
  }

  Vector4D& Vector4D::operator*=(const Vector4D& rhs)
  {
    this->x_ *= rhs.x_;
    this->y_ *= rhs.y_;
    this->z_ *= rhs.z_;
    this->w_ *= rhs.w_;
    return *this;
  }

  Vector4D& Vector4D::operator/=(const Vector4D& rhs)
  {
    this->x_ /= rhs.x_;
    this->y_ /= rhs.y_;
    this->z_ /= rhs.z_;
    this->w_ /= rhs.w_;
    return *this;
  }

  void Vector4D::axpy(const double mult, const Vector4D& rhs)
  {
    this->x_ += mult * rhs.x_;
    this->y_ += mult * rhs.y_;
    this->z_ += mult * rhs.z_;
    this->w_ += mult * rhs.w_;
  }

  double Vector4D::dot_product_with(const Vector4D& rhs) const
  {
    return this->x_ * rhs.x_ + this->y_ * rhs.y_ + this->z_ * rhs.z_ + this->w_ * rhs.w_;
  }

  void Vector4D::multiply(const Vector4D& rhs, Vector4D& lhs)
  {
    lhs.x_ = x_ * rhs.x_;
    lhs.y_ = y_ * rhs.y_;
    lhs.z_ = z_ * rhs.z_;
    lhs.w_ = w_ * rhs.w_;
  }

  void Vector4D::print(std::ostream & os) const {
    os << x_ << ", " << y_ << ", " << z_ << ", " << w_ << std::endl;
  }
}  // namespace test
