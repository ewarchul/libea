#pragma once

#include <concepts>

namespace libea::common {
template <typename T>
concept Numeric = std::integral<T> || std::floating_point<T>;
}
