#include <libcmaes/cmaes.h>
#include <libcmaes/esoptimizer.h>

#include <iostream>
#include <libea/solvers/cmaes/cmaes.hpp>
#include <libea/termination/predefined.hpp>
#include <libea/termination/termination.hpp>

using namespace libcmaes;

namespace chr = std::chrono;

using namespace libea::termination;


template <typename Parameters, typename Solutions>
const auto custom_stops =
    termination_criteria<Parameters, Solutions>{.criteria_ = {
                                                    max_iter<Parameters, Solutions>,
                                                }};



FitFunc fsphere = [](const double* x, const int N) {
  double val = 0.0;
  for (int i = 0; i < N; i++) val += x[i] * x[i];
  return val;
};

constexpr auto sphere_fn = [](const auto& x) {
  double result{0.0};
  for (std::size_t i{0}; i < x.size(); ++i) {
    result += x[i] * x[i];
  }
  return result;
};

auto run_libcmaes(const std::size_t N) {
  const auto dim{N};
  std::vector<double> x0(dim, 95);
  double sigma = 0.1;
  CMAParameters<> cmaparams(x0, sigma);
  std::cout << "lambda = " << cmaparams.lambda() << std::endl;
  cmaparams.set_algo(CMAES_DEFAULT);
  CMASolutions cmasols = cmaes<>(fsphere, cmaparams);
  std::cout << "best solution: " << cmasols << std::endl;
}

auto run_eacxx(const auto N) {
  namespace cma_es = libea::solvers::cmaes;
  const auto dim{N};
  Eigen::VectorXd x0 = Eigen::VectorXd::Ones(dim) * 95;
  auto params = cma_es::parameters{dim};
  auto stops = custom_stops<cma_es::parameters, cma_es::solutions>;
  params.max_iter_ = 1400;
  fmt::print("{}", params);
  auto csa_cma_es = cma_es::csa_cmaes(x0, sphere_fn, params, stops);
  auto start = chr::system_clock::now();
  const auto ret = csa_cma_es.solve();
  auto end = chr::system_clock::now();
  fmt::print("CSA-CMA-ES optimization results: {}\n", ret);
}

auto main() -> int {
  {
    fmt::print("Running libcmaes\n");
    auto start = chr::system_clock::now();
   // run_libcmaes(100);
    auto end = chr::system_clock::now();
    fmt::print("LIBCMAES -> {}ms\n", chr::duration_cast<chr::milliseconds>(end - start).count());
  }
  {
    fmt::print("Running eacxx\n");
    auto start = chr::system_clock::now();
    run_eacxx(100u);
    auto end = chr::system_clock::now();
    fmt::print("LIBCMAES -> {}ms\n", chr::duration_cast<chr::milliseconds>(end - start).count());
  }

  return 0;
}
