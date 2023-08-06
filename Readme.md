# libeacpp [work-in-progress]

**{libeacpp}** is [**{blaze}**-ingly](https://bitbucket.org/blaze-lib/blaze/src/master/) fast header-only library providing evolutionary optimizition
solvers written in C++20.

## Requirements

## Build & install

## Solvers

## Example 


```cpp
#include <libea/solvers/cmaes/cmaes.hpp>
#include <libea/termination/predefined.hpp>

namespace cma_es = libea::solvers::cmaes;

// evaluation function
constexpr auto sphere_fn = [](const auto& x) {
  double result{0.0};
  for (int i{0}; i < x.size(); ++i) {
    result += x[i] * x[i];
  }
  return result;
};

auto main() -> int {
  constexpr int dim{10}; // number of decision variables
  constexpr blaze::StaticVector<double, dim> x0(95); // init point
  const auto params = cma_es::parameters{dim};  // create CMA-ES parameters

  const auto stops =  // predefined termination conditions
    libea::termination::predefined_termination_criteria<cma_es::parameters, cma_es::solutions>; 

  const auto csa_cma_es = // vanila CMA-ES, e.s, Cumulative Step-size Adaptation (CSA) Covariance Matrix Adaptation Evolution Strategy 
    cma_es::csa_cmaes(x0, sphere_fn, params, stops); 

  const auto solutions = csa_cma_es.run(); // optimization results

  return 0;
}
```


## Benchmarks


