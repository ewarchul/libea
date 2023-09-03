#pragma once

#include <fmt/core.h>
#include <libea/common/consts.h>
#include <libea/common/types.h>

#include <Eigen/Core>
#include <iostream>

namespace libea::solvers::cmaes {

struct parameters {
  parameters() = default;

  parameters(types::u32_t dimension) : dim_{dimension} {
    max_iter_ = types::as<types::u64_t>(std::floor(max_fevals_ / lambda_));
    tolx_ = consts::DEFAULT_TOLX_MUL * sigma_;

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
  types::u64_t max_iter_{};
  types::u32_t lambda_{4 * dim_};
  types::u32_t mu_{types::as<types::u32_t>(std::floor(types::as<double>(lambda_) / 2))};

  double tolx_{};
  double tol_fitness_{consts::DEFAULT_FITNESS_TOL};
  double target_fitness_{std::numeric_limits<double>::min()};

  double tol_cov_cond_{};

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

template <> struct fmt::formatter<libea::solvers::cmaes::parameters> {
  constexpr auto parse(format_parse_context& ctx) -> format_parse_context::iterator {
    auto it = ctx.begin(), end = ctx.end();
    if (it != end && *it != '}') fmt::detail::throw_format_error("invalid format");
    return it;
  }

  auto format(const libea::solvers::cmaes::parameters& p, fmt::format_context& ctx) -> fmt::format_context::iterator {
    return fmt::format_to(ctx.out(),
                          "\n--\ndim = {};\nn.max.fevals = {};\nn.max.iter = {};\nxtol "
                          "= {};\ninit.pop.size = {};\ninit.sigma = {};\n--\n",
                          p.dim_, p.max_fevals_, p.max_iter_, p.tolx_, p.lambda_, p.sigma_);
  }
};
