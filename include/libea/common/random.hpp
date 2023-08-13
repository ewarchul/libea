#pragma once

#include <EigenRand/EigenRand>
#include <random>

namespace {
std::random_device rd{};
Eigen::Rand::P8_mt19937_64 urng{rd()};
}  // namespace

namespace libea::common {

// @TODO change type aliases
inline auto rnorm(const auto rows, const auto cols)  {
  return Eigen::Rand::normal<Eigen::MatrixXd>(rows, cols, urng);
}
}  // namespace libea::common
