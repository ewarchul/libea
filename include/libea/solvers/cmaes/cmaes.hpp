#pragma once

#include <libea/solvers/cmaes/solver.hpp>
#include <libea/updaters/sigma/csa.hpp>
#include <libea/updaters/sigma/idle.hpp>

namespace libea::solvers::cmaes {
constexpr auto idle_cmaes = [](const auto& x0, auto&& fn, const auto& params, const auto& terminations) {
  auto idle_updater = libea::updaters::sigma::idle{};
  return solver{x0, fn, idle_updater, params, terminations};
};

constexpr auto csa_cmaes = [](const auto& x0, auto&& fn, const auto& params, const auto& terminations) {
  auto idle_updater = libea::updaters::sigma::csa{params};
  return solver{x0, fn, idle_updater, params, terminations};
};


}
