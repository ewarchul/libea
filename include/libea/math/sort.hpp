#pragma once

#include <range/v3/algorithm/sort.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/take.hpp>

namespace libea::math {
// @TODO remove repetitions
auto trunc_sort_ind(auto &&fitness_vals, const auto trunc_num) {
  auto zipped = fitness_vals | ranges::views::enumerate | ranges::to_vector;
  ranges::sort(zipped, [](const auto &r1, const auto &r2) {
    return r1.second < r2.second;
  });
  auto indices = zipped | ranges::views::take(trunc_num) |
                 ranges::views::keys | ranges::to<std::vector<int>>;
  auto fitness = zipped | ranges::views::take(trunc_num) |
                 ranges::views::values | ranges::to<std::vector<int>>;
  return std::make_pair(indices, fitness);
}
} // namespace libea::math
