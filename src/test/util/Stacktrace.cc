/*
 * (C) Copyright 2021-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
/*
 * This test deliberately produces a stacktrace.
 */

#include "oops/util/Stacktrace.h"

#include <iostream>

int main(int argc, char **argv) {
  std::cout << "Stacktrace:\n" << util::stacktrace_current() << std::endl;
  return 0;
}

