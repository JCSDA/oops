/*
 * (C) Crown Copyright 2025-2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string_view>
#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "oops/util/Printable.h"

namespace oops {
namespace mpi {

/// \brief Utility class to infer information regarding the MPI color
///        grouping on a given split communicator.
///
/// \details For example, given a split communicator where there are 6 ranks,
///          and two color groups 0 and 1 split as follows:
///
///          {0, 0, 1, 1, 1, 1}
///
///          Say we are on MPI color 1. This class will return:
///
///          ```
///          colorInfo.color() => 1
///          colorInfo.roots() => {0, 2}
///          colorInfo.numColors() => 2
///          ```
///
///          **Assumptions**
///          1) Color groups are defined in order. i.e, a color group {1, 2, 0} is problematic
///             as there is no way of telling that PE 0 is color 1 without that prior knowledge.
///             What ColorInfo would produce in that case is: {0, 1, 2}.
///
class ColorInfo : public util::Printable {
 public:
  /// \brief Constructor.
  /// \param subComm
  ///    A split MPI subcommunicator.
  /// \param parentComm
  ///    The MPI communicator that was split to produce the subComm.
  ColorInfo(const eckit::mpi::Comm& subComm,
            const eckit::mpi::Comm& parentComm = eckit::mpi::comm("world"));

  /// \brief Constructor.
  /// \param subCommName
  ///    A split MPI subcommunicator name.
  /// \param parentCommName
  ///    The name of the MPI communicator that was split to produce the subComm.
  ColorInfo(const std::string_view subCommName,
            const std::string_view parentCommName = "world");

  /// \brief The color this rank resides in in the split communicator.
  size_t color() const { return color_; }

  /// \brief The parent ranks of the root rank of each color.
  /// \details For example, with an MPI communicator split across 6 ranks: {0, 0, 1, 1, 1, 1},
  ///          the roots are {0, 2}.
  const std::vector<size_t>& roots() const { return parentRoots_; }

  /// \brief The number of color groups within this sub communicator.
  size_t numColors() const { return parentRoots_.size(); }

  void print(std::ostream &) const override;

 private:
  static const std::string classname() { return "ColorInfo"; }

  size_t color_;

  /// \brief Vector containing the root rank of each color, in terms of the parent rank numbering.
  std::vector<size_t> parentRoots_;

  const std::string parentName_;
  const std::string subName_;
};

}  // namespace mpi
}  // namespace oops
