#pragma once

#include <blaze/Math.h>

namespace libea::updaters::sigma {

// Cumulative Step-size adaptation rule 
// https://ieeexplore.ieee.org/document/542381
struct csa {
  static constexpr bool requires_fitness_function{false};
  csa(const auto& params) : chin_{params.chin_}, csigma_{params.csigma_} {
    damping_factor_ = 1 + 2 * blaze::max(0, blaze::sqrt((params.mueff_ - 1) / (params.dim_ + 1)) - 1) + csigma_;
  }

  auto operator()(double sigma, const auto& solutions) -> double {
    return sigma * blaze::exp(blaze::norm(solutions.psigma_) / chin_ - 1) * csigma_ / damping_factor_;
  }

  double chin_;
  double csigma_;
  double damping_factor_{};
};

}  // namespace libea::updaters::sigma
