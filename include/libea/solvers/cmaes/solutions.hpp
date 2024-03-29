#pragma once

#include <libea/common/types.h>
#include <libea/math/generators.hpp>
#include <limits>
#include <vector>

namespace libea::solvers::cmaes {

struct solutions {
  solutions() = default;

  solutions(types::u32_t dimension) : dim_{dimension} {
    psigma_ = blaze::zero<types::number_t>(dim_);
    pcov_ = blaze::zero<types::number_t>(dim_);

    eigen_vecs = math::diag(dim_, 1);
    BD = eigen_vecs * math::diag(dim_, 1);
    cov_mat = BD * blaze::trans(BD);
  }

  auto inc_feval(const types::u32_t lambda) { fevals_ += lambda; }
  auto inc_iter() { ++iter_; }

  auto update_best(const auto& current_best, const auto fitness) {
    if (fitness < best_so_far_fitness_) {
      best_so_far_ = current_best;
      best_so_far_fitness_ = fitness;
    }
    best_log_.push_back(current_best);
    best_log_fitness.push_back(fitness);
  }

  types::u32_t dim_{};
  unsigned iter_{};
  unsigned fevals_{};

  double sigma_{1};

  types::dvec_t mean_{};
  types::dmat_t population_{};

  types::dvec_t best_so_far_{};
  types::number_t best_so_far_fitness_{std::numeric_limits<types::number_t>{}.max()};
  std::vector<types::dvec_t> best_log_{};
  std::vector<types::number_t> best_log_fitness{};

  bool hsig_{true};
  types::dvec_t pcov_;
  types::dvec_t psigma_;
  types::dmat_t eigen_vecs;
  types::dmat_t BD;
  blaze::SymmetricMatrix<types::dmat_t> cov_mat{};
};

}  // namespace libea::solvers::cmaes
