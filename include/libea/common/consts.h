#pragma once

#include <libea/common/types.h>

namespace libea::consts {
static constexpr unsigned MAX_ITER_MUL{10000};
static constexpr unsigned MAX_FEVALS_MUL{10000};
static constexpr unsigned LAMBDA_MUL{4};
static constexpr types::number_t DEFAULT_TOL{10e-6};
static constexpr types::number_t DEFAULT_LOWER_BOUND{-100.0};
static constexpr types::number_t DEFAULT_UPPER_BOUND{100.0};
static constexpr types::number_t DEFAULT_STEP_SIZE{1};
} // namespace libea::consts
