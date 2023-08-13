#pragma once

#include <Eigen/Core>
#include <blaze/math/DynamicVector.h>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/take.hpp>
#include <libea/common/types.h>

#include <iostream>

namespace libea::math {

auto evaluate(const auto &population, auto &&fitness_fn) {
  auto colview = population.colwise();
  auto result = ranges::views::transform(colview, fitness_fn) | ranges::to_vector;
  Eigen::VectorXd fitness_values = Eigen::Map<Eigen::VectorXd>(result.data(), result.size(), 1);
  return fitness_values;

}

} // namespace libea::math
