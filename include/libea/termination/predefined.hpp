#pragma once

#include <libea/common/types.h>

#include <cmath>
#include <iostream>
#include <libea/termination/termination.hpp>

namespace libea::termination {

template <typename P, typename S> using criterium = termination_criteria<P, S>::value_type;

template <typename P, typename S>
const auto max_fevals = criterium<P, S>{"[MaxFevals] maximal number of function evaluations reached",
                                        [](const auto& p, const auto& s) { return s.fevals_ > p.max_fevals_; }};

template <typename P, typename S>
const auto max_iter = criterium<P, S>{"[MaxIter] maximal number of iteration reached",
                                      [](const auto& p, const auto& s) { return s.iter_ > p.max_iter_; }};

template <typename P, typename S>
const auto tol_x =
    criterium<P, S>{"[TolX] standard deviation below tolerance in all coordinates", [](const auto& p, const auto& s) {
                      const auto std_cond = (s.eigen_values_.array() < p.tolx_).sum() == p.dim_;
                      const auto pcov_cond = ((s.sigma_ * s.pcov_).array() < p.tolx_).sum() == p.dim_;
                      return std_cond and pcov_cond;
                    }};

template <typename P, typename S>
const auto tol_fitness = criterium<P, S>{"[TolFitness] optimal function value approximated", [](const auto& p, const auto& s) {
                                           return std::abs(s.best_so_far_fitness_ - p.target_fitness_) < p.tol_fitness_;
                                         }};

template <typename P, typename S>
const auto condition_number = criterium<P, S>{"[ConditionCov] condition number of covariance matrix exceeds assumed threshold",
                                              [](const auto& p, const auto& s) {
                                                if (s.iter_ < 1) {
                                                  return false;
                                                }
                                                const auto cond_num = s.eigen_values_.maxCoeff() / s.eigen_values_.minCoeff();
                                                return cond_num < p.tol_cov_cond_;
                                              }};

template <typename P, typename S>
const auto indef_cov_mat = criterium<P, S>{"[IndefCovMat] covariance matrix is not numerically positive definite",
                                           [](const auto&, const auto& s) { return s.indef_cov_mat_; }};

template <typename P, typename S>
const auto no_axis_effect = criterium<P, S>{
    "[NoEffectAxis] Addition of 10% of step size multiplier does not change mean value", [](const auto& p, const auto& s) {
      if (s.iter_ < 1) {
        return false;
      }
      const auto index = (s.iter_ % p.dim_) + 1;
      const auto eigen_vec = s.eigen_vecs_.col(index);
      const auto eigen_val = s.eigen_values_[index];
      const auto diff_vec = s.mean_ - (s.mean_ + 0.1 * s.sigma_ * std::sqrt(eigen_val) * eigen_vec);
      return diff_vec.norm() < std::numeric_limits<double>::min();
    }};

template <typename P, typename S>
const auto no_coord_effect =
    criterium<P, S>{"[NoEffectCoord] Addition of 20% of step size multiplier in any coordinate does not change mean value",
                    [](const auto&, const auto& s) {
                      if (s.iter_ < 1) {
                        return false;
                      }
                      const types::dvec_t shifted_mean = s.mean_.array() + 0.2 * s.sigma_;
                      const auto diff = s.mean_ - shifted_mean;
                      return diff.norm() < std::numeric_limits<double>::min();
                    }};

template <typename Parameters, typename Solutions>
const auto predefined_termination_criteria =
    termination_criteria<Parameters, Solutions>{.criteria_ = {
                                                    max_fevals<Parameters, Solutions>,
                                                    tol_x<Parameters, Solutions>,
                                                    tol_fitness<Parameters, Solutions>,
                                                    condition_number<Parameters, Solutions>,
                                                    indef_cov_mat<Parameters, Solutions>,
                                                    no_axis_effect<Parameters, Solutions>,
                                                    no_coord_effect<Parameters, Solutions>,
                                                }};

}  // namespace libea::termination
