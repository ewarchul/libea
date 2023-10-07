#pragma once

#include <Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h>
#include <libea/Math.h>

#include <libea/common/types.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <EigenRand/EigenRand>
#include <concepts>
#include <exception>
#include <iostream>
#include <libea/solvers/cmaes/parameters.hpp>
#include <libea/solvers/cmaes/solutions.hpp>
#include <libea/termination/termination.hpp>
#include <utility>

#include "libea/common/random.hpp"

namespace libea::solvers::cmaes {

template <typename SigmaUpdater, typename FitnessFunction> class solver {
 public:
  solver(const types::dvec_t& x0, SigmaUpdater sigma_updater, FitnessFunction fitness_fn, parameters params,
         termination::termination_criteria<parameters, solutions> stops)
      : sigma_updater_{std::move(sigma_updater)},
        fitness_fn_{std::move(fitness_fn)},
        params_{std::move(params)},
        stops_{std::move(stops)} {
    solutions_ = solutions{types::as<types::u32_t>(x0.size()), params.lambda_};
    solutions_.mean_ = x0;
  }

  [[nodiscard]] auto solve() -> solutions {
    while (!terminate()) {
      auto [pop, diffs] = ask();
      const auto fn_vals = math::evaluate(pop, fitness_fn_);
      solutions_.inc_feval(params_.lambda_);
      solutions_.fitness_values_ = fn_vals;

      solutions_.prev_population_ = std::move(solutions_.population_);
      solutions_.population_ = std::move(pop);

      auto [selected_indices, selected_fn_vals] = math::trunc_sort_ind(fn_vals, params_.mu_);
      auto selected_pop = solutions_.population_(Eigen::all,
                                                 selected_indices);  // blaze::columns(pop, selected_indices);
      auto selected_diffs = diffs(Eigen::all, selected_indices);

      auto current_best = selected_pop(Eigen::all, 0);
      solutions_.update_best(current_best, selected_fn_vals[0]);

      tell(std::move(selected_pop), std::move(selected_diffs));

      solutions_.inc_iter();
    }

    return solutions_;
  }

 private:
  [[nodiscard]] auto ask() {
    types::dmat_t diffs = common::rnorm(params_.dim_, params_.lambda_);
    std::cout << "DIFFS" << std::endl;
    std::cout << diffs << std::endl;
    std::cout << "DIFFS END\n";
    types::dmat_t population = solutions_.mean_.replicate(1, params_.lambda_) + solutions_.sigma_ * solutions_.BD_mat_ * diffs;

    return std::make_pair(population, diffs);
  }

  auto tell(auto&& selected_pop, auto&& selected_diffs) -> void {
    solutions_.mean_ = selected_pop * params_.weights_;

    const auto zmean = selected_diffs * params_.weights_;

    solutions_.psigma_ =
        params_.psigma_decay_factor_ * solutions_.psigma_ + params_.psigma_coeff_ * (solutions_.eigen_vecs_ * zmean);

    const auto hsig_exp = 2 * solutions_.fevals_ / params_.lambda_;
    solutions_.hsig_ = (solutions_.psigma_.norm() / std::sqrt(1.0 - std::pow(params_.psigma_decay_factor_, hsig_exp)) /
                        params_.chin_) < params_.hsig_coeff_;

    solutions_.pcov_ =
        params_.pcov_decay_factor_ * solutions_.pcov_ + solutions_.hsig_ * params_.pcov_coeff_ * (solutions_.BD_mat_ * zmean);

    update_cma(std::move(selected_diffs));
    update_sigma();

    eigen_solver_.compute(solutions_.cov_mat_);
    if (eigen_solver_.info() != Eigen::Success) [[unlikely]] {
      solutions_.indef_cov_mat_ = true;
    } else {
      solutions_.eigen_values_ = eigen_solver_.eigenvalues().real().unaryExpr([](const auto& x) { return std::sqrt(x); });
      solutions_.eigen_vecs_ = eigen_solver_.eigenvectors().real();
      solutions_.BD_mat_ = solutions_.eigen_vecs_ * solutions_.eigen_values_.asDiagonal();
    }
  }

  auto update_sigma() -> void {
    if constexpr (SigmaUpdater::requires_fitness_function) {
      solutions_.sigma_ = sigma_updater_(solutions_.sigma_, solutions_, fitness_fn_);
    } else {
      solutions_.sigma_ = sigma_updater_(solutions_.sigma_, solutions_);
    }
  }

  auto update_cma(auto&& selz) -> void {
    const auto BDz = solutions_.BD_mat_ * selz;

    const auto empirical_cov = (1.0 - params_.ccov_) * solutions_.cov_mat_;
    const auto rank_mu_cov = params_.ccov_ * (1.0 / params_.cmu_) *
                             (solutions_.pcov_ * solutions_.pcov_.transpose() +
                              (1.0 - solutions_.hsig_) * params_.cc_ * (2.0 - params_.cc_) * solutions_.cov_mat_);
    const auto rank_one_cov = params_.ccov_ * (1 - 1 / params_.cmu_) * BDz * params_.weights_diag_ * BDz.transpose();

    solutions_.cov_mat_ = empirical_cov + rank_mu_cov + rank_one_cov;
  }

  [[nodiscard]] auto terminate() -> bool {
    const auto [terminate, msg] = stops_.check(params_, solutions_);
    if (terminate) {
      solutions_.stop_msg_ = msg;
    }
    return terminate;
  }

  SigmaUpdater sigma_updater_;
  FitnessFunction fitness_fn_;
  parameters params_;
  solutions solutions_;
  termination::termination_criteria<parameters, solutions> stops_;
  Eigen::SelfAdjointEigenSolver<types::dmat_t> eigen_solver_;
};

}  // namespace libea::solvers::cmaes
