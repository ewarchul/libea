#pragma once

#include <libea/solvers/cmaes/solver.hpp>
#include <libea/updaters/sigma/csa.hpp>
#include <libea/updaters/sigma/idle.hpp>
#include <libea/updaters/sigma/ppmf.hpp>

namespace libea::solvers::cmaes {
constexpr auto idle_cmaes = [](const auto& x0, auto&& fn, const auto& params,
                               const auto& terminations) {
  auto sigma_updater = libea::updaters::sigma::idle{};
  return solver{x0, sigma_updater, fn, params, terminations};
};

constexpr auto csa_cmaes = [](const auto& x0, auto&& fn, const auto& params,
                              const auto& terminations) {
  auto sigma_updater = libea::updaters::sigma::csa{params};
  return solver{x0, sigma_updater, fn, params, terminations};
};

constexpr auto ppmf_cmaes = [](const auto& x0, auto&& fn, const auto& params,
                               const auto& terminations) {
  auto sigma_updater = libea::updaters::sigma::ppmf{2.0, 0.1};
  return solver{x0, sigma_updater, fn, params, terminations};
};

}  // namespace libea::solvers::cmaes
