#pragma once

#include <random>

namespace {
std::random_device rd{};
std::mt19937 gen{rd()};
} // namespace

namespace libea::common {

auto rnorm() -> double {
  auto iso_normal_dist = std::normal_distribution<>{0, 1};
  return iso_normal_dist(gen);
}
} // namespace libea::common
