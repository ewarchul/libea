#pragma once

#include <libea/common/types.h>

#include <Eigen/Core>
#include <iostream>
#include <limits>
#include <vector>

namespace libea::solvers::cmaes {

struct solutions {
  solutions() = default;

  solutions(types::u32_t dimension, types::u32_t lambda) : dim_{dimension}, lambda_{lambda} {
    psigma_ = types::dvec_t::Zero(dim_);
    pcov_ = types::dvec_t::Zero(dim_);

    population_ = types::dmat_t::Zero(dim_, lambda_);
    prev_population_ = types::dmat_t::Zero(dim_, lambda_);

    auto e_vec = types::dvec_t::Ones(dim_);
    eigen_vecs = e_vec.asDiagonal();

    types::dmat_t f = e_vec.asDiagonal();
    BD = f * e_vec.asDiagonal();

    cov_mat = BD * BD.transpose();
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
  types::u64_t iter_{};
  types::u64_t dummy_iter_{};
  types::u64_t fevals_{};
  types::u32_t lambda_{};
  double sigma_{1};

  types::dvec_t mean_{};

  types::dmat_t population_{};
  types::dmat_t prev_population_{};

  types::dvec_t fitness_values_{};
  types::dvec_t best_so_far_{};
  double best_so_far_fitness_{std::numeric_limits<double>{}.max()};
  std::vector<types::dvec_t> best_log_{};
  std::vector<double> best_log_fitness{};

  bool hsig_{true};
  types::dvec_t pcov_;
  types::dvec_t psigma_{};
  types::dmat_t eigen_vecs{};
  types::dmat_t BD{};
  types::dmat_t cov_mat{};
};

}  // namespace libea::solvers::cmaes
