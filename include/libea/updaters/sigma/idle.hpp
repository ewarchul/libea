#pragma once

namespace libea::updaters::sigma {

// Idle step-size adaptation rule for testing purposes
struct idle {
  static constexpr bool requires_fitness_function{false};
  auto operator()(double step_size, const auto& solutions) -> double { return step_size * param_; }
  double param_{1};
};

}  // namespace libea::updaters::sigma
