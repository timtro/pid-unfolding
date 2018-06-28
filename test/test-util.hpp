#pragma once

#include <type_traits>
#include <vector>

namespace util {

  template <typename T, typename A, typename F>
  auto fmap(F f, const std::vector<T, A> &v) {
    using MapedToType = std::invoke_result_t<F, decltype(v.at(0))>;
    std::vector<MapedToType> result;
    result.reserve(v.size());

    for (const auto &each : v) result.push_back(f(each));

    return result;
  }

}  // namespace util
