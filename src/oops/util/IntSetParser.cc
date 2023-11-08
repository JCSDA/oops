/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/IntSetParser.h"

#include <algorithm>

#include <sstream>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Function to split string on delimiter

std::vector<std::string> splitString(const std::string& str, char delim)
{
  std::vector<std::string> result;
  std::istringstream s(str);
  std::string substr;
  while (std::getline(s, substr, delim)) {
    result.push_back(substr);
  }
  return result;
}

// -----------------------------------------------------------------------------

/// Function to convert string to channel number. String to integer conversion in standard
//  routines (eg, std::stoi) will sometimes throw when non-numeric characers are
//  in the input string instead of throwing an invalid_argument exception. This routine
//  checks for non-numeric characters and returns a -1 (channel numbers should be positive)
//  if such characters exist. The blank (' ') is included in the find_first_not_of below
//  since the splitString routine (above) is deliminting only on a comma (',') which leaves
//  blanks in its result.

int stringToChanNum(const std::string& str) {
  // Abort if input string contains non-valid characters
  if (str.find_first_not_of("0123456789 ") != std::string::npos) {
    throw eckit::BadValue("Input string contains non-numeric characters");
  }

  int chnum;
  std::istringstream ss(str);
  ss >> chnum;
  return chnum;
}

// -----------------------------------------------------------------------------

/// Function to parse integers (supports commas for separating integers
//  and integers ranges and dashes for integer ranges).
//  For example: 1-5, 9, 13-45
//  Returns a std::set, no need to sort or remove duplicates and find/insert are in log(n)
//  Supports -1 as [-1] and not [0-1]; other negative integers are not parsed

std::set<int> parseIntSet(const std::string & str) {
  std::set<int> set_ints;

// split string by commas to get individual integers or ranges
  std::vector<std::string> ranges = splitString(str, ',');

  if (ranges == std::vector<std::string>(1, "-1")) {
    set_ints.insert(-1);
  } else {
    for (std::size_t irange = 0; irange < ranges.size(); irange++) {
      // split the element by dashes (in case it is a range)
      std::vector<std::string> range = splitString(ranges[irange], '-');
      ASSERT((range.size() == 1) || (range.size() == 2));
      // add a single channel
      if (range.size() == 1) {
        // add a single channel
        set_ints.insert(stringToChanNum(range[0]));
      } else if (range.size() == 2) {
        // add a range
        int start = stringToChanNum(range[0]);
        int stop  = stringToChanNum(range[1]);
        for (int ch = start; ch <= stop; ch++) {
          set_ints.insert(ch);
        }
      }
    }
  }
  return set_ints;
}

// -----------------------------------------------------------------------------
}  // namespace oops
