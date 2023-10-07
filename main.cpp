#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/ranges.h>

#include <Eigen/Core>
#include <chrono>
#include <iostream>
#include <libea/solvers/cmaes/cmaes.hpp>
#include <libea/termination/predefined.hpp>
#include <range/v3/algorithm/for_each.hpp>

namespace cma_es = libea::solvers::cmaes;
namespace chr = std::chrono;

constexpr auto sphere_fn = [](const auto& x) {
  double result{0.0};
  for (std::size_t i{0}; i < x.size(); ++i) {
    result += x[i] * x[i];
  }
  return result;
};

auto save_solutions(const auto& filename, const auto& solutions) {
  auto output = fmt::output_file(filename);
  output.print("t, best_so_far, best_current, sigma, dim\n");
  const auto iter = solutions.iter_;
  ranges::for_each(ranges::views::iota(0u, iter), [&](const auto t) {
    output.print("{},{:0.8f},{:0.8f}, {:0.15f},{}\n", t, solutions.best_so_far_fitness_history[t], solutions.best_fitness_history[t], solutions.sigma_history[t], solutions.dim_);
  });
}

auto main() -> int {
  constexpr int dim{40};
  Eigen::VectorXd x0 = Eigen::VectorXd::Ones(dim) * 95;
  auto params = cma_es::parameters{dim};
  auto stops = libea::termination::predefined_termination_criteria<cma_es::parameters, cma_es::solutions>;
  params.max_iter_ = 1400;
  auto csa_cma_es = cma_es::csa_cmaes(x0, sphere_fn, params, stops);
  auto ppmf_cma_es = cma_es::ppmf_cmaes(x0, sphere_fn, params, stops);

  // {
  //   fmt::print("Running CSA-CMA-ES with following parameters: {}\n", params);
  //   auto start = chr::system_clock::now();
  //   const auto ret = csa_cma_es.solve();
  //   auto end = chr::system_clock::now();
  //
  //   fmt::print("elapsed time -> {}ms\n", chr::duration_cast<chr::milliseconds>(end - start).count());
  //   fmt::print("CSA-CMA-ES optimization results: {}\n", ret);
  //
  //   save_solutions("csa-cma-es-spehere-100D.csv", ret);
  // }

{
    fmt::print("Running PPMF-CMA-ES with following parameters: {}\n", params);
    auto start = chr::system_clock::now();
    const auto ret = ppmf_cma_es.solve();
    auto end = chr::system_clock::now();

    fmt::print("elapsed time -> {}ms\n", chr::duration_cast<chr::milliseconds>(end - start).count());
    fmt::print("PPMF-CMA-ES optimization results: {}\n", ret);

    save_solutions("ppmf-cma-es-spehere-100D.csv", ret);
  }


  return 0;
}
