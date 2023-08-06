#pragma once

#include <blaze/Math.h>

#include <iostream>
#include <range/v3/algorithm/count_if.hpp>
#include <stdexcept>

namespace libea::updaters::sigma {

// Previous Population Midpoint Fitness step-size adaptation rule
// https://ieeexplore.ieee.org/abstract/document/9504829
struct ppmf {
  static constexpr bool requires_fitness_function{true};
  ppmf(double damping_factor, double target_probability) : damp_factor{damping_factor}, target_prob{target_probability} {
    if (damping_factor < 0.0) {
      throw std::invalid_argument("Damping factor cannot be negative!");
    }
    if (target_probability < 0.0) {
      throw std::invalid_argument("Target probability cannot be negative!");
    }
  }
  auto operator()(double sigma, auto& solutions, auto&& fitness_func) -> double {
    const auto prev_midpoint = blaze::mean<blaze::columnwise>(solutions.prev_population_);
    const auto prev_midpoint_fitness = fitness_func(prev_midpoint);
    solutions.inc_feval(1);
    double success_prob = ranges::count_if(solutions.fitness_values_, [&](auto fn_val) { return fn_val < prev_midpoint_fitness; }) / solutions.lambda_;
    std::cout << success_prob << std::endl;
    if (solutions.iter_ == 10) {
      exit(1);
    }
    return sigma * blaze::exp(damp_factor * blaze::abs(success_prob - target_prob) / (1 - target_prob));
  }

  double damp_factor{};
  double target_prob{};
};

}  // namespace libea::updaters::sigma
