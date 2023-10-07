#pragma once

#include <blaze/Math.h>
#include <fmt/core.h>
#include <libea/common/types.h>

#include <iostream>
#include <range/v3/algorithm/count_if.hpp>
#include <stdexcept>

namespace libea::updaters::sigma {

// Previous Population Midpoint Fitness step-size adaptation rule
// https://ieeexplore.ieee.org/abstract/document/9504829

// @TODO investigate why it stops to work with bigger dimensions
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
    std::cout << "NEWW!!\n";
    const auto prev_midpoint = solutions.prev_population_.rowwise().mean();
    const auto prev_midpoint_fitness = fitness_func(prev_midpoint);
    solutions.inc_feval(1);
    const double success_prob =
        ranges::count_if(solutions.fitness_values_, [&](auto fn_val) {
          return fn_val < prev_midpoint_fitness;
        }) / solutions.lambda_;

    std::cout << "++++++\n";
    std::cout << solutions.population_ << std::endl;
    std::cout << "++++++\n";


       // std::cout << "########\n";
    // std::cout << solutions.fitness_values_ << std::endl;
    // std::cout << "fffffffffffff\n";
    // std::cout << prev_midpoint << std::endl;
    // std::cout << "===\n";
    // fmt::print("q(x_avg_prev) = {:5.10f}\n", prev_midpoint_fitness);
    // std::cout << "+++\n";
    // fmt::print("P_success = {:5.10f}\n", success_prob);
    // std::cout << "dddd\n";
    // std::cout << solutions.lambda_ << std::endl;


    auto new_sigma =  sigma * std::exp(damp_factor * (success_prob - target_prob) / (1 - target_prob));
    // fmt::print("old sigma  =  {}\n", sigma);
    // fmt::print("damp-factor = {}\n", damp_factor);
    // fmt::print("P_target = {}\n", target_prob);
    // fmt::print("new sigma  =  {}\n", new_sigma);
    if (solutions.iter_ == 5) { 
      exit(1);
    }
    return new_sigma;
  }

  double damp_factor{};
  double target_prob{};
};

}  // namespace libea::updaters::sigma
