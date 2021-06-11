/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/wildcard.h"

namespace util {

// An std::string-based implementation of the algorithm from Robert van Engelen,
// "Fast String Matching with Wildcards, Globs, and Gitignore-Style Globs - How Not to Blow it Up",
// https://www.codeproject.com/Articles/5163931/Fast-String-Matching-with-Wildcards-Globs-and-Giti
bool matchesWildcardPattern(const std::string &string, const std::string &pattern) {
  const std::size_t stringSize = string.size();
  const std::size_t patternSize = pattern.size();

  if (patternSize == 0)
    return stringSize == 0;

  std::size_t stringPos = 0;
  std::size_t patternPos = 0;
  std::size_t stringBacktrackPos = std::string::npos;
  std::size_t patternBacktrackPos = std::string::npos;

  while (stringPos != stringSize) {
    if (patternPos != patternSize && pattern[patternPos] == '*') {
      // Found a new * wildcard. Set new positions to use in case of backtracking
      stringBacktrackPos = stringPos;
      patternBacktrackPos = ++patternPos;
    } else if (patternPos != patternSize &&
               (pattern[patternPos] == string[stringPos] || pattern[patternPos] == '?')) {
      // Matched a single character
      ++stringPos;
      ++patternPos;
    } else {
      // Mismatch or end of pattern
      if (patternBacktrackPos == std::string::npos) {
        // There have been no * wildcards -- we can't backtrack.
        // The string doesn't match the pattern.
        return false;
      }
      // Backtrack: try to assign one more character to the most recent * wildcard.
      stringPos = ++stringBacktrackPos;
      patternPos = patternBacktrackPos;
    }
  }

  // Ignore any trailing * wildcards.
  while (patternPos != patternSize && pattern[patternPos] == '*')
    patternPos++;

  // Have we matched the complete pattern?
  return patternPos == patternSize;
}

}  // namespace util
