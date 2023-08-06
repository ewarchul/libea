#pragma once

#include <blaze/math/dense/DynamicMatrix.h>
#include <blaze/math/dense/DynamicVector.h>

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
using number_t = double;

using dvec_t = blaze::DynamicVector<number_t>;
using dmat_t = blaze::DynamicMatrix<number_t>;
} // namespace libea::types
