#include <iostream>

#include <catch/catch.hpp>

#include "../lib/pid.hpp"
#include "test-util.hpp"

using Real_t = double;
using RealSignalVector = std::vector<SignalPt<Real_t>>;
using CtrlState = PIDState<Real_t>;

const auto now = chrono::steady_clock::now();

constexpr auto scanl =
    [](const auto f, const auto x0,
       const RealSignalVector &data) -> std::vector<CtrlState> {
  auto accumulator = x0;
  std::vector<decltype(accumulator)> xs;
  xs.reserve(data.size());

  for (const auto &datum : data) {
    accumulator = f(datum, accumulator);
    xs.push_back(accumulator);
  };

  return xs;
};

TEST_CASE(
    "A proportional ctrler with kp = 1, given a constant error signal with "
    "value 0 should produce a constant control signal of value 0 and "
    "accumulate no error.") {
  const auto foldable_p00 = pid_algebra<double>(1., 0., 0.);

  RealSignalVector ones{{now + 1s, 0}, {now + 2s, 0}, {now + 3s, 0}};

  CtrlState init{now, 0., 0., 0.};

  auto result = scanl(foldable_p00, init, ones);

  REQUIRE(util::fmap(result, [](const CtrlState &s) { return s.time; })
          == std::vector<std::decay_t<decltype(now)>>{now + 1s, now + 2s,
                                                      now + 3s});
  REQUIRE(util::fmap(result, [](const CtrlState &s) { return s.ctrlVal; })
          == std::vector<double>{0., 0., 0.});
  REQUIRE(util::fmap(result, [](const CtrlState &s) { return s.error; })
          == std::vector<double>{0., 0., 0.});
  REQUIRE(util::fmap(result, [](const CtrlState &s) { return s.errSum; })
          == std::vector<double>{0., 0., 0.});
}

TEST_CASE(
    "An integral controller with ki = 1, given a constant error signal of "
    "value 1, should accumulate error with a proportionally growing ctrlVal.") {
  const auto foldable_0i0 = pid_algebra<double>(0., 1., 0.);

  RealSignalVector ones{{now + 1s, 1}, {now + 2s, 1}, {now + 3s, 1}};

  CtrlState init{now, 0., 0., 0.};

  auto result = scanl(foldable_0i0, init, ones);

  REQUIRE(util::fmap(result, [](const CtrlState &s) { return s.time; })
          == std::vector<std::decay_t<decltype(now)>>{now + 1s, now + 2s,
                                                      now + 3s});
  REQUIRE(util::fmap(result, [](const CtrlState &s) { return s.ctrlVal; })
          == std::vector<double>{1., 2., 3.});
  REQUIRE(util::fmap(result, [](const CtrlState &s) { return s.error; })
          == std::vector<double>{1., 1., 1.});
  REQUIRE(util::fmap(result, [](const CtrlState &s) { return s.errSum; })
          == std::vector<double>{1., 2., 3.});
}

TEST_CASE(
    "A derivative controller with kd = 1, given a unit amplitude sawtooth "
    "should demonstrate the appropriate derivative/difference output. The "
    "error accumulation should also reflect the sawtooth signal.") {
  const auto foldable_00d = pid_algebra<double>(0., 0., 1.);

  RealSignalVector ones{{now + 1s, 0}, {now + 2s, 1}, {now + 3s, 0}};

  CtrlState init{now, 0., 0., 0.};

  auto result = scanl(foldable_00d, init, ones);

  REQUIRE(util::fmap(result, [](const CtrlState &s) { return s.time; })
          == std::vector<std::decay_t<decltype(now)>>{now + 1s, now + 2s,
                                                      now + 3s});
  REQUIRE(util::fmap(result, [](const CtrlState &s) { return s.ctrlVal; })
          == std::vector<double>{0., 1., -1.});
  REQUIRE(util::fmap(result, [](const CtrlState &s) { return s.error; })
          == std::vector<double>{0., 1., 0.});
  REQUIRE(util::fmap(result, [](const CtrlState &s) { return s.errSum; })
          == std::vector<double>{0., 1., 1.});
}
