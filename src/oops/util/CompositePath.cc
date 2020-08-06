/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/CompositePath.h"

namespace util {

CompositePath::CompositePath(char separator)
  : separator_(separator), path_(1, separator)
{}

void CompositePath::push(const std::string &component) {
  if (path_.size() > 1) {
    // We're not at the root, so we need to append a separator to the previous component
    path_ += separator_;
  }
  path_ += component;
}

void CompositePath::pop(size_t componentLength) {
  // Remove the last component
  path_.resize(path_.size() - componentLength);
  if (path_.size() > 1) {
    // We're not at the root, so remove also the final separator
    path_.resize(path_.size() - 1);
  }
}

PathComponent::PathComponent(CompositePath &path, const std::string &newComponent)
  : path_(path), length_(newComponent.size())
{
  path_.push(newComponent);
}

PathComponent::~PathComponent() {
  path_.pop(length_);
}

}  // namespace util
