#include <fmt/core.h>

#include <chrono>
#include <iostream>
#include <libea/solvers/cmaes/cmaes.hpp>
#include <libea/termination/predefined.hpp>

namespace cma_es = libea::solvers::cmaes;
namespace chr = std::chrono;

constexpr auto sphere_fn = [](const auto& x) {
  double result{0.0};
  for (int i{0}; i < x.size(); ++i) {
    result += x[i] * x[i];
  }
  return result;
};

auto main() -> int {
  constexpr int dim{10};
  blaze::StaticVector<double, dim> x0(95);
  auto params = cma_es::parameters{dim};
  auto stops = libea::termination::predefined_termination_criteria<cma_es::parameters, cma_es::solutions>;

  auto ppmf_cma_es = cma_es::ppmf_cmaes(x0, sphere_fn, params, stops);

  auto start = chr::system_clock::now();
  const auto ret = ppmf_cma_es.solve();
  auto end = chr::system_clock::now();

  fmt::print("elapsed time -> {}ms\n", chr::duration_cast<chr::milliseconds>(end - start).count());
  std::cout << ret.best_so_far_fitness_ << std::endl;

  return 0;
}