#pragma once

#include <libea/common/consts.h>
#include <libea/common/types.h>

#include <Eigen/Core>
#include <iostream>

namespace libea::solvers::cmaes {

struct parameters {
  parameters() = default;

  parameters(types::u32_t dimension) : dim_{dimension} {
    weights_ = std::log(mu_ + 1) - types::dvec_t::LinSpaced(mu_, 1, mu_).array().log();
    weights_ = weights_ / weights_.sum();
    weights_diag_ = weights_.asDiagonal();

    mueff_ = std::pow(weights_.sum(), 2) / weights_.unaryExpr([](const auto x) { return x * x; }).sum();
    cmu_ = mueff_;
    ccov_ = (1.0 / cmu_) * 2.0 / std::pow(dim_ + 1.4, 2) +
            (1.0 - 1.0 / cmu_) * ((2.0 * cmu_ - 1.0) / (std::pow((dim_ + 2), 2.0) + 2.0 * cmu_));

    csigma_ = (mueff_ + 2) / (dim_ + mueff_ + 3);
    chin_ = std::sqrt(dim_) * (1.0 - (1.0 / (4.0 * dim_)) + (1.0 / (21.0 * dim_ * dim_)));
    psigma_coeff_ = std::sqrt(csigma_ * (2.0 - csigma_) * mueff_);
    psigma_decay_factor_ = 1.0 - csigma_;
    pcov_decay_factor_ = 1.0 - cc_;
    pcov_coeff_ = std::sqrt(cc_ * (2.0 - cc_) * mueff_);
  }

  types::u32_t dim_;
  double hsig_coeff_{(1.4 + 2.0 / (dim_ + 1.0))};

  types::u64_t max_fevals_{consts::MAX_FEVALS_MUL * dim_};
  types::u32_t lambda_{consts::LAMBDA_MUL * dim_};
  types::u32_t mu_{types::as<types::u32_t>(std::floor(types::as<double>(lambda_) / 2))};

  double xtol_{consts::DEFAULT_TOL};
  double lower_bound_{consts::DEFAULT_LOWER_BOUND};
  double upper_bound_{consts::DEFAULT_UPPER_BOUND};

  double sigma_{consts::DEFAULT_STEP_SIZE};
  double cc_{4.0 / (dim_ + 4.0 + 0.0)};
  double csigma_{};

  double psigma_decay_factor_{};
  double psigma_coeff_{};

  double pcov_decay_factor_{};
  double pcov_coeff_{};

  double ccov_{};
  double chin_{};
  double cmu_{};
  double mueff_{};

  types::dvec_t weights_{};
  types::dmat_t weights_diag_{};
};
}  // namespace libea::solvers::cmaes
