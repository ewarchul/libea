#pragma once

#include <libea/common/types.h>

namespace libea::consts {
static constexpr unsigned MAX_FEVALS_MUL{1000};
static constexpr unsigned LAMBDA_MUL{4};
static constexpr double DEFAULT_TOLX_MUL{10e-12};
static constexpr double DEFAULT_FITNESS_TOL{10e-12};
static constexpr double DEFAULT_LOWER_BOUND{-100.0};
static constexpr double DEFAULT_UPPER_BOUND{100.0};
static constexpr double DEFAULT_STEP_SIZE{1};
}  // namespace libea::consts
