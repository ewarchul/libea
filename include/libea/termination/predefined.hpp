#pragma once

#include <libea/termination/termination.hpp>

namespace libea::termination {
template <typename Parameters, typename Solutions>
const auto predefined_termination_criteria = termination_criteria<Parameters, Solutions>{

    .criteria_ = {{"max_fevals", [](const auto& p, const auto& s) { return s.fevals_ > p.max_fevals_; }},
                  {"max_iter", [](const auto& p, const auto& s) { return s.fevals_ > p.max_fevals_; }},
                  {"stop_fitness", [](const auto& p, const auto& s) { return s.fevals_ > p.max_fevals_; }}}};
}
