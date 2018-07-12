#pragma once

#include <chrono>
#include <iostream>
#include <type_traits>
#include <vector>
#include <optional>
#include <utility>

#include <sodium/sodium.h>
#include <catch/catch.hpp>

namespace chrono = std::chrono;

namespace util {

  template <typename T, typename A, typename F>
  auto fmap(F f, const std::vector<T, A> &v) {
    using MapedToType = std::invoke_result_t<F, decltype(v.at(0))>;
    std::vector<MapedToType> result;
    result.reserve(v.size());

    for (const auto &each : v) result.push_back(f(each));

    return result;
  }

  // unfold : ( (A → optional<pair<A, B>>), A ) → vector<B>
  template <typename F, typename A>
  auto unfold(F f, A a0) {
    // will fail if `f` doesn't return a pair when given an A.
    using B = decltype(std::declval<std::invoke_result_t<F, A>>()->second);
    static_assert(std::is_same_v<std::invoke_result_t<F, A>,
                                 std::optional<std::pair<A, B>>>);

    std::vector<B> bs;

    auto resultAB = f(a0);
    if (resultAB)
      bs.push_back(resultAB->second);
    else
      return bs;

    while (1) {
      resultAB = f(resultAB->first);
      if (resultAB)
        bs.push_back(resultAB->second);
      else
        return bs;
    }
  }

  template <typename Clock = chrono::steady_clock>
  inline constexpr auto double_to_duration(double t) {
    return chrono::duration_cast<typename Clock::duration>(
        chrono::duration<double>(t));
  }

  template <typename A>
  inline constexpr auto unchrono_sec(A t) {
    chrono::duration<double> tAsDouble = t;
    return tAsDouble.count();
  }

  template <typename T>
  bool compareVectors(const std::vector<T> &a, const std::vector<T> &b,
                      T margin) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++) {
      if (a[i] != Approx(b[i]).margin(margin)) {
        std::cout << a[i] << " Should == " << b[i] << std::endl;
        return false;
      }
    }
    return true;
  }
}  // namespace util
