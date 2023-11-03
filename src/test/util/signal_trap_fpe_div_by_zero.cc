/*
 * (C) Copyright 2021-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
/*
 * This test deliberately produces a floating point exception. Not using eckit
 * because our fpe handlers always produce an abort() call. ctest will expect
 * this test to always fail.
 */

#include "oops/util/signal_trap.h"

#include <iostream>

int main(int argc, char **argv) {
  util::trap_sigfpe(1);
  std::cout << "Test to trigger an fpe division by zero." << std::endl;

  for (int i=10; i >= 0; i--) {
    std::cout << "10 / " << i << " = ";
    double res = 10. / static_cast<double>(i);
    std::cout << res << std::endl;
  }
  return 0;
}

