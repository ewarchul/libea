#pragma once

#include <Eigen/Core>
#include <Eigen/src/Core/Matrix.h>
#include <blaze/math/DynamicVector.h>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/take.hpp>
#include <fmt/core.h>
#include <libea/common/types.h>

#include <iostream>

namespace libea::math {

// @TODO rewrite in idomatic Eigen, i.e., with unaryExp
auto evaluate(const auto &population, auto &&fitness_fn) -> Eigen::VectorXd {
  auto colview = population.colwise();
  auto result = ranges::views::transform(colview, fitness_fn) | ranges::to<std::vector<double>>;
  Eigen::VectorXd fitness_values = Eigen::Map<Eigen::VectorXd>(result.data(), result.size(), 1);
  return fitness_values;
}

} // namespace libea::math
