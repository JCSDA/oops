/*
 * (C) Crown copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_ALGORITHMS_H_
#define OOPS_UTIL_ALGORITHMS_H_

#include <type_traits>
#include <vector>

namespace util {

/// \brief Returns the vector containing the results of calling the function `op` on each element
/// of the vector `inputs`.
template <typename Input, typename UnaryOperation>
std::vector<typename std::invoke_result<UnaryOperation, const Input &>::type> transformVector(
    const std::vector<Input> &inputs, const UnaryOperation &op)
{
  typedef typename std::invoke_result<UnaryOperation, const Input&>::type Output;
  std::vector<Output> outputs;
  outputs.reserve(inputs.size());
  for (const Input &input : inputs)
    outputs.push_back(op(input));
  return outputs;
}

}  // namespace util

#endif  // OOPS_UTIL_ALGORITHMS_H_
