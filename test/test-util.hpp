#pragma once

#include <vector>

namespace util {

  template <typename T, typename A, typename F>
  auto get_each(std::vector<T, A> v, F accessor) {
    using AccessedType = std::invoke_result_t<F, decltype(v.at(0))>;

    std::vector<AccessedType> result;
    result.reserve(v.size());

    for (const auto &each : v) result.push_back(accessor(each));

    return result;
  }
}  // namespace util
