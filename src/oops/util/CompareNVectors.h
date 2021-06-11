/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_COMPARENVECTORS_H_
#define OOPS_UTIL_COMPARENVECTORS_H_

#include <stddef.h>
#include <vector>

namespace oops {

  /// \brief Tests whether at least one vector is empty.
  /// This is the base case for one vector.
  template <typename T>
    bool anyVectorEmpty(const std::vector <T> &vec)
  {
    return vec.empty();
  }

  /// \brief Tests whether at least one vector is empty.
  /// This is the recursive case that accepts an arbitrary number of vectors
  /// using a variadic template.
  template <typename T, typename... Args>
    bool anyVectorEmpty(const std::vector <T> &vec1,
                        const Args&... vecs)
  {
    return vec1.empty() || anyVectorEmpty(vecs...);
  }

  /// \brief Tests whether all vectors have the same size.
  /// This is the trivial base case for one vector. The result always evaluates to true.
  template <typename T>
    bool allVectorsSameSize(const std::vector <T> &vec)
  {
    return true;
  }

  /// \brief Tests whether all vectors have the same size.
  /// This is the base case for two vectors.
  template <typename T1, typename T2>
    bool allVectorsSameSize(const std::vector <T1> &vec1,
                            const std::vector <T2> &vec2)
  {
    return vec1.size() == vec2.size();
  }

  /// \brief Tests whether all vectors have the same size.
  /// This is the recursive case that accepts an arbitrary number of vectors
  /// using a variadic template.
  /// The first vector in the list is compared against all subsequent ones.
  template <typename T1, typename T2, typename... Args>
    bool allVectorsSameSize(const std::vector <T1> &vec1,
                            const std::vector <T2> &vec2,
                            const Args&... vecs)
  {
    return vec1.size() == vec2.size() && allVectorsSameSize(vec1, vecs...);
  }

  /// \brief Tests whether all vectors have an expected size (N).
  /// This is the base case for one vector.
  template <typename T>
    bool allVectorsExpectedSize(const size_t &N,
                                const std::vector <T> &vec)
  {
    return vec.size() == N;
  }

  /// \brief Tests whether all vectors have an expected size (N).
  /// This is the recursive case that accepts an arbitrary number of vectors
  /// using a variadic template.
  template <typename T, typename... Args>
    bool allVectorsExpectedSize(const size_t &N,
                                const std::vector <T> &vec,
                                const Args&... vecs)
  {
    return vec.size() == N && allVectorsExpectedSize(N, vecs...);
  }

  /// \brief Tests whether all filled (non-empty) vectors have an expected size (N).
  /// This is the base case for one vector.
  template <typename T>
    bool allNonEmptyVectorsExpectedSize(const size_t &N,
                                        const std::vector <T> &vec)
  {
    return vec.empty() || vec.size() == N;
  }

  /// \brief Tests whether all filled (non-empty) vectors have an expected size (N).
  /// This is the recursive case that accepts an arbitrary number of vectors
  /// using a variadic template.
  template <typename T, typename... Args>
    bool allNonEmptyVectorsExpectedSize(const size_t &N,
                                        const std::vector <T> &vec,
                                        const Args&... vecs)
  {
    return (vec.empty() || vec.size() == N) && allNonEmptyVectorsExpectedSize(N, vecs...);
  }

  /// \brief Tests whether all vectors have the same nonzero size.
  /// This is the base case for one vector.
  template <typename T>
    bool allVectorsSameNonZeroSize(const std::vector <T> &vec)
  {
    return !vec.empty();
  }

  /// \brief Tests whether all vectors have the same nonzero size.
  /// This is the base case for two vectors.
  template <typename T1, typename T2>
    bool allVectorsSameNonZeroSize(const std::vector <T1> &vec1,
                                   const std::vector <T2> &vec2)
  {
    return !vec1.empty() && vec1.size() == vec2.size();
  }

  /// \brief Tests whether all vectors have the same nonzero size.
  /// This is the recursive case that accepts an arbitrary number of vectors
  /// using a variadic template.
  /// The first vector in the list is compared against all subsequent ones.
  template <typename T1, typename T2, typename... Args>
    bool allVectorsSameNonZeroSize(const std::vector <T1> &vec1,
                                   const std::vector <T2> &vec2,
                                   const Args&... vecs)
  {
    return !vec1.empty() && vec1.size() == vec2.size() && allVectorsSameSize(vec1, vecs...);
  }
}  // namespace oops
#endif  // OOPS_UTIL_COMPARENVECTORS_H_
