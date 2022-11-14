#pragma once

#include <vector>

#include <sodium/sodium.h>

namespace util {

  template <typename T>
  auto make_listener(const sodium::stream<T> &s) {
    auto sRecord = std::make_shared<std::vector<T>>();
    const auto s_unlisten =
        s.listen([sRecord](auto x) { sRecord->push_back(x); });
    return std::make_pair(s_unlisten, sRecord);
  }

}  // namespace util
