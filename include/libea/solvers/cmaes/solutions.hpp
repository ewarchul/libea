#pragma once

#include <fmt/core.h>
#include <libea/common/types.h>

#include <Eigen/Core>
#include <iostream>
#include <limits>
#include <vector>

namespace libea::solvers::cmaes {

struct solutions {
  solutions() = default;

  solutions(types::u32_t dimension, types::u32_t lambda)
      : dim_{dimension}, lambda_{lambda} {
    psigma_ = types::dvec_t::Zero(dim_);
    pcov_ = types::dvec_t::Zero(dim_);

    population_ = types::dmat_t::Zero(dim_, lambda_);
    population_.fill(10000000);

    prev_population_ = types::dmat_t::Zero(dim_, lambda_);
    prev_population_.fill(1000000);



    auto e_vec = types::dvec_t::Ones(dim_);
    eigen_vecs_ = e_vec.asDiagonal();

    types::dmat_t f = e_vec.asDiagonal();
    BD_mat_ = f * e_vec.asDiagonal();

    cov_mat_ = BD_mat_ * BD_mat_.transpose();
  }

  auto inc_feval(const types::u32_t lambda) { fevals_ += lambda; }
  auto inc_iter() { ++iter_; }

  auto update_best(const auto& current_best, const auto fitness) {
    if (fitness < best_so_far_fitness_) {
      best_so_far_ = current_best;
      best_so_far_fitness_ = fitness;
    }
    sigma_history.push_back(sigma_);
    best_param_history.push_back(current_best);
    best_fitness_history.push_back(fitness);
    best_so_far_fitness_history.push_back(best_so_far_fitness_);
  }

  types::u32_t dim_{};
  types::u64_t iter_{};
  types::u64_t fevals_{};
  std::string stop_msg_{};

  types::u32_t lambda_{};
  double sigma_{1};

  types::dvec_t mean_{};

  types::dmat_t population_{};
  types::dmat_t prev_population_{};

  types::dvec_t fitness_values_{};
  types::dvec_t best_so_far_{};
  double best_so_far_fitness_{std::numeric_limits<double>{}.max()};
  std::vector<types::dvec_t> best_param_history{};
  std::vector<double> best_fitness_history{};
  std::vector<double> best_so_far_fitness_history{};
  std::vector<double> sigma_history{};

  bool hsig_{true};
  types::dvec_t pcov_;
  types::dvec_t psigma_{};
  types::dmat_t eigen_vecs_{};
  types::dvec_t eigen_values_{};
  types::dmat_t eigen_vals_diag_{};
  types::dmat_t BD_mat_{}; // eigen_vecs * eigen_vals_diag_
  types::dmat_t cov_mat_{};
  bool indef_cov_mat_{false};
};

}  // namespace libea::solvers::cmaes

template <> struct fmt::formatter<libea::solvers::cmaes::solutions> {
  constexpr auto parse(format_parse_context& ctx)
      -> format_parse_context::iterator {
    auto it = ctx.begin(), end = ctx.end();
    if (it != end && *it != '}')
      fmt::detail::throw_format_error("invalid format");
    return it;
  }

  auto format(const libea::solvers::cmaes::solutions& s,
              fmt::format_context& ctx) -> fmt::format_context::iterator {
    return fmt::format_to(ctx.out(),
                          "\n--\nn.iter = {};\nn.fevals = {};\nbest_fit = "
                          "{};\nstop.msg = {};\n--\n",
                          s.iter_, s.fevals_, s.best_so_far_fitness_,
                          s.stop_msg_);
  }
};
