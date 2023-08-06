#pragma once

#include <blaze/math/DynamicVector.h>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/take.hpp>

namespace libea::math {

auto evaluate(const auto &population, auto &&fitness_fn) {
  blaze::DynamicVector<double> fitness_values(
      population.columns());
  for (auto column :
       ranges::views::iota(0) | ranges::views::take(population.columns())) {
    const auto &col_view = blaze::column(population, column);
    fitness_values[column] = fitness_fn(col_view);
  }
  return fitness_values;
}
} // namespace libea::math
