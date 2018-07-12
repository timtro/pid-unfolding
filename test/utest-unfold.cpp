#include <type_traits>
#include <utility>

#include <catch/catch.hpp>

#include "../include/util/util.hpp"

TEST_CASE(
    "Anamorphic fold should return an empty list if the coalgebra returns a "
    "nullopt given the seed.") {
  auto coalg = [](int) -> std::optional<std::pair<int, int>> { return {}; };

  REQUIRE(util::unfold(coalg, 41) == std::vector<int>{});
}

TEST_CASE(
    "Given a coalg : int → optional<pair<int, bool>>, and an int seed, the "
    "unfold function should return a vector<bool>.") {
  auto coalg = [](int) -> std::optional<std::pair<int, bool>> { return {}; };

  REQUIRE(util::unfold(coalg, 41) == std::vector<bool>{});
}

TEST_CASE(
    "Anamorphic unfold should count down from seed to zero (inclusive) given "
    "this coalgebra") {
  auto coalg = [](int current) -> std::optional<std::pair<int, int>> {
    auto oneSmaller = current - 1;
    if (current < 0)
      return {};
    else
      return {{oneSmaller, current}};
  };

  REQUIRE(util::unfold(coalg, 10)
          == std::vector{10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0});
}
