#include <fmt/core.h>

#include <Eigen/Core>
#include <chrono>
#include <iostream>
#include <libea/solvers/cmaes/cmaes.hpp>
#include <libea/termination/predefined.hpp>

namespace cma_es = libea::solvers::cmaes;
namespace chr = std::chrono;

constexpr auto sphere_fn = [](const auto& x) {
  double result{0.0};
  for (std::size_t i{0}; i < x.size(); ++i) {
    result += x[i] * x[i];
  }
  return result;
};

auto main() -> int {
  constexpr int dim{100};
  Eigen::VectorXd x0 = Eigen::VectorXd::Ones(dim) * 95;
  auto params = cma_es::parameters{dim};
  auto stops =
      libea::termination::predefined_termination_criteria<cma_es::parameters,
                                                          cma_es::solutions>;
  auto csa_cma_es = cma_es::csa_cmaes(x0, sphere_fn, params, stops);

  {
    fmt::print("Running CSA-CMA-ES with following parameters: {}\n", params);
    auto start = chr::system_clock::now();
    const auto ret = csa_cma_es.solve();
    auto end = chr::system_clock::now();

    fmt::print("elapsed time -> {}ms\n",
               chr::duration_cast<chr::milliseconds>(end - start).count());
    fmt::print("CSA-CMA-ES optimization results: {}\n", ret);
  }

  return 0;
}
