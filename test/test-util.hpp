#pragma once

#include <vector>
#include <type_traits>

namespace util {

  template <typename T, typename A, typename F>
  auto fmap(const std::vector<T, A>& v, F f) {
    using MapedToType = std::invoke_result_t<F, decltype(v.at(0))>;
    std::vector<MapedToType> result;
    result.reserve(v.size());

    for (const auto &each : v) result.push_back(f(each));

    return result;
  }

}// namespace util
