#pragma once

namespace libea::updaters::sigma {

struct idle {
  auto operator()(double step_size, const auto& solutions) -> double { return step_size * param_; }

  double param_{1};
};

} // namespace libea::updaters::sigma
