#pragma once


#include <range/v3/range/conversion.hpp>
#include <range/v3/view/iota.hpp>

#include <libea/common/consts.h>
#include <libea/common/types.h>
#include <libea/math/generators.hpp>

namespace libea::solvers::cmaes {

struct parameters {
  parameters() = default;

  parameters(types::u32_t dimension) : dim_{dimension} {
    weights_ = blaze::log(mu_ + 1) - blaze::log(blaze::linspace(mu_, types::as<types::u32_t>(1), mu_));
    weights_ = weights_ / blaze::sum(weights_);
    weights_diag_ = math::diag(weights_);

    mueff_ = blaze::pow(blaze::sum(weights_), 2) / blaze::sum(blaze::pow(weights_, 2));
    cmu_ = mueff_;
    ccov_ = (1.0 / cmu_) * 2.0 / blaze::pow(dim_ + 1.4, 2) +
            (1.0 - 1.0 / cmu_) * ((2.0 * cmu_ - 1.0) / (blaze::pow((dim_ + 2), 2.0) + 2.0 * cmu_));

    csigma_ = (mueff_ + 2) / (dim_ + mueff_ + 3);
    chin_ = blaze::sqrt(dim_) * (1.0 - (1.0 / (4.0 * dim_)) + (1.0 / (21.0 * dim_ * dim_)));
    psigma_coeff_ = blaze::sqrt(csigma_ * (2.0 - csigma_) * mueff_);
    psigma_decay_factor_ = 1.0 - csigma_;
    pcov_decay_factor_ = 1.0 - cc_;
    pcov_coeff_ = blaze::sqrt(cc_ * (2.0 - cc_) * mueff_);
  }

  types::u32_t dim_;
  types::number_t hsig_coeff_{(1.4 + 2.0 / (dim_ + 1.0))};

  types::u64_t max_fevals_{consts::MAX_FEVALS_MUL * dim_};
  types::u32_t lambda_{consts::LAMBDA_MUL * dim_};
  types::u32_t mu_{types::as<types::u32_t>(blaze::floor(types::as<types::number_t>(lambda_) / 2))};

  types::number_t xtol_{consts::DEFAULT_TOL};
  types::number_t lower_bound_{consts::DEFAULT_LOWER_BOUND};
  types::number_t upper_bound_{consts::DEFAULT_UPPER_BOUND};

  types::number_t sigma_{consts::DEFAULT_STEP_SIZE};
  types::number_t cc_{4.0 / (dim_ + 4.0 + 0.0)};
  types::number_t csigma_{};

  types::number_t psigma_decay_factor_{};
  types::number_t psigma_coeff_{};

  types::number_t pcov_decay_factor_{};
  types::number_t pcov_coeff_{};

  types::number_t ccov_{};
  types::number_t chin_{};
  types::number_t cmu_{};
  types::number_t mueff_{};
  types::dvec_t weights_{};
  types::dmat_t weights_diag_{};
};
}  // namespace ew_cmaes
