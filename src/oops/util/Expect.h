/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_EXPECT_H_
#define OOPS_UTIL_EXPECT_H_

#include <sstream>
#include <string>

// Recent versions of the JCSDA fork of eckit define the EXPECT_EQUAL macro in eckit/testing/Test.h.
// This definition is repeated below for compatibility with older versions of eckit.
#ifndef EXPECT_EQUAL

// IMPORTANT: To use the macro below, it is also necessary to include "eckit/testing/Test.h",
// after defining ECKIT_TESTING_SELF_REGISTER_CASES if needed.

// Provides more informative output on failure than the raw EXPECT() macro.
#define EXPECT_EQUAL(expr, expected) \
    do { \
        if (!((expr) == (expected))) { \
            std::stringstream str; \
            str << ("EXPECT condition '" #expr " == " #expected "' failed. ") \
                << "(Received: " << expr << "; expected: " << expected << ")"; \
            throw eckit::testing::TestException(str.str(), Here()); \
        } \
    } while (false)

#endif  // EXPECT_EQUAL

/// Check if \p expr throws an exception with a message containing the string \p msg.
#define EXPECT_THROWS_MSG(expr, msg)                                                         \
  do {                                                                                       \
    bool exceptionWithMsgThrown = false;                                                     \
    try {                                                                                    \
      expr;                                                                                  \
    }                                                                                        \
    catch (const std::exception &ex) {                                                       \
      if (std::strstr(ex.what(), msg) != nullptr)                                            \
        exceptionWithMsgThrown = true;                                                       \
    }                                                                                        \
    if (!exceptionWithMsgThrown)                                                             \
      throw eckit::testing::TestException("Expected exception with message '" +              \
                                          std::string(msg) +                                 \
                                          "' not thrown in: " #expr, Here());                \
  } while (false)

#endif  // OOPS_UTIL_EXPECT_H_
