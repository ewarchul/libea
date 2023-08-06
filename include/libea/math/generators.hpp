#pragma once

#include <libea/common/concepts.h>
#include <libea/common/random.hpp>
#include <range/v3/range/concepts.hpp>
#include <blaze/Math.h>


namespace libea::math {
inline auto diag(auto size, common::Numeric auto scalar) {
  return blaze::generate(size, size, [scalar](auto row, auto col) { return row == col ? scalar : 0; });
}

inline auto diag(const ranges::range auto& vector) {
  return blaze::generate(vector.size(), vector.size(), [&vector](auto row, auto col) { return row == col ? vector[row] : 0; });
}

inline auto diag(const ranges::range auto& vector, auto&& op) {
  return blaze::generate(vector.size(), vector.size(), [&](auto row, auto col) { return row == col ? op(vector[row]) : 0; });
}

inline auto random(auto rows, auto cols) {
  return blaze::generate(rows, cols, [](auto r, auto c) { return common::rnorm(); });
}

}  // namespace math

