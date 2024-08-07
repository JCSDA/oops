/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "test/assimilation/Vector2D.h"

namespace test {

  Vector2D::Vector2D(const double& x,
                     const double& y)
    : x_(x), y_(y)
  {}

  Vector2D::Vector2D(const Vector2D& rhs,
                     const bool copy)
    : x_(rhs.x_), y_(rhs.y_)
  {
    // NB the argument 'copy' is present in order to match the signature
    // required by various minimiser functions.
  }

  Vector2D& Vector2D::operator=(const Vector2D& rhs)
  {
    this->x_ = rhs.x_;
    this->y_ = rhs.y_;
    return *this;
  }

  Vector2D& Vector2D::operator+=(const Vector2D& rhs)
  {
    this->x_ += rhs.x_;
    this->y_ += rhs.y_;
    return *this;
  }

  Vector2D& Vector2D::operator-=(const Vector2D& rhs)
  {
    this->x_ -= rhs.x_;
    this->y_ -= rhs.y_;
    return *this;
  }

  Vector2D& Vector2D::operator*=(const double mult)
  {
    this->x_ *= mult;
    this->y_ *= mult;
    return *this;
  }

  Vector2D& Vector2D::operator*=(const Vector2D& rhs)
  {
    this->x_ *= rhs.x_;
    this->y_ *= rhs.y_;
    return *this;
  }

  Vector2D& Vector2D::operator/=(const Vector2D& rhs)
  {
    this->x_ /= rhs.x_;
    this->y_ /= rhs.y_;
    return *this;
  }

  void Vector2D::axpy(const double mult, const Vector2D& rhs)
  {
    this->x_ += mult * rhs.x_;
    this->y_ += mult * rhs.y_;
  }

  double Vector2D::dot_product_with(const Vector2D& rhs) const
  {
    return this->x_ * rhs.x_ + this->y_ * rhs.y_;
  }

  void Vector2D::multiply(const Vector2D& rhs, Vector2D& lhs)
  {
    lhs.x_ = x_ * rhs.x_;
    lhs.y_ = y_ * rhs.y_;
  }

  void Vector2D::print(std::ostream & os) const {
    os << x_ << ", " << y_ << std::endl;
  }
}  // namespace test
