#pragma once

#include <fmt/core.h>
#include <functional>
#include <string>
#include <unordered_map>

namespace libea::termination {

using termination_name = std::string;

struct termination_status {
  bool terminate_{false};
  std::string msg_;
};

template <typename Parameters, typename Solutions> struct termination_criteria {
  using termination_callback = std::function<bool(const Parameters&, const Solutions&)>;

  auto check(const auto& params, const auto& sols) -> termination_status {
    for (auto&& [name, terminate_cb] : criteria_) {
      if (terminate_cb(params, sols)) {
        return termination_status{.terminate_ = true,
                                  .msg_ = fmt::format("Solver is terminated with the following "
                                                      "termination criteria: {}",
                                                      name)};
      }
    }
    return {};
  }

  auto add_critiera(const termination_name& name, const termination_callback& cb) -> void { criteria_.try_emplace(name, cb); }

  auto remove_critiera(const termination_name& name) {
    if (criteria_.contains(name)) {
      criteria_.erase(name);
    }
  }

  std::unordered_map<termination_name, termination_callback> criteria_;
};

}  // namespace libea::termination
