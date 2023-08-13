#pragma once

#include <Eigen/Core>

namespace libea::types {

constexpr auto to_underlying(auto input)
    -> std::underlying_type_t<decltype(input)> {
  return static_cast<std::underlying_type_t<decltype(input)>>(input);
}

template <typename OutputType> constexpr auto as(auto input) -> OutputType {
  return static_cast<OutputType>(input);
}

using u8_t = uint8_t;
using u16_t = uint16_t;
using u32_t = uint32_t;
using u64_t = uint64_t;

using dvec_t = Eigen::VectorXd;
using dmat_t = Eigen::MatrixXd;
} // namespace libea::types
