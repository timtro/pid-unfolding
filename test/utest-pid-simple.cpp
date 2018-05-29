#include <catch/catch.hpp>
#include <iostream>
#include "../lib/pid.hpp"

using Real_t           = double;
using RealSignalVector = std::vector<SignalPt<Real_t>>;
using CtrlState        = PIDState<Real_t>;

const auto now    = chrono::steady_clock::now();
const auto oneSec = 1s;

template <typename T, typename A, typename F>
auto get_each(std::vector<T, A> v, F accessor) {
  using AccessedType = std::invoke_result_t<F, decltype(v.at(0))>;
  std::vector<AccessedType> result;
  result.reserve(v.size());
  for (const auto &each : v) {
    result.push_back(accessor(each));
    std::vector<T, A> result(v.size());
  }
  return result;
}

constexpr auto list_fold =
    [](const auto f, const auto x0,
       const RealSignalVector &data) -> std::vector<CtrlState> {
  auto accumulator = x0;
  std::vector<decltype(accumulator)> xs;
  xs.reserve(data.size());

  for (const auto &datum : data) {
    accumulator = f(accumulator, datum);
    xs.push_back(accumulator);
  };
  return xs;
};

TEST_CASE(
    "A proportional ctrler with kp = 1, given a constant signal with value 1 "
    "and a setpoint of 1 should produce a constant control signal of value 0, "
    "accumulate no error. It also shouldn't report error or change its "
    "setpoint.") {
  const auto foldable_p00 = pid_algebra<double>(1.f, 0.f, 0.f);

  RealSignalVector ones{{now + 1s, 1}, {now + 2s, 1}, {now + 3s, 1}};

  CtrlState init{now, 1.f, 0.f, 0.f, 0.f};

  auto result = list_fold(foldable_p00, init, ones);

  REQUIRE(get_each(result, [](const CtrlState &s) { return s.time; })
          == std::vector<std::decay_t<decltype(now)>>{now + 1s, now + 2s,
                                                      now + 3s});
  REQUIRE(get_each(result, [](const CtrlState &s) { return s.ctrlVal; })
          == std::vector<double>{0.f, 0.f, 0.f});
  REQUIRE(get_each(result, [](const CtrlState &s) { return s.error; })
          == std::vector<double>{0.f, 0.f, 0.f});
  REQUIRE(get_each(result, [](const CtrlState &s) { return s.errSum; })
          == std::vector<double>{0.f, 0.f, 0.f});
  REQUIRE(get_each(result, [](const CtrlState &s) { return s.setPt; })
          == std::vector<double>{1.f, 1.f, 1.f});
}

TEST_CASE(
    "An integral controller with ki = 1, given a constant signal of value 1 "
    "and a setpoint of zero should accumulate error with a proportionally "
    "growing ctrlVal.") {
  const auto foldable_0i0 = pid_algebra<double>(0.f, 1.f, 0.f);

  RealSignalVector ones{{now + 1s, 1}, {now + 2s, 1}, {now + 3s, 1}};

  CtrlState init{now, 0.f, 0.f, 0.f, 0.f};

  auto result = list_fold(foldable_0i0, init, ones);

  REQUIRE(get_each(result, [](const CtrlState &s) { return s.time; })
          == std::vector<std::decay_t<decltype(now)>>{now + 1s, now + 2s,
                                                      now + 3s});
  REQUIRE(get_each(result, [](const CtrlState &s) { return s.ctrlVal; })
          == std::vector<double>{1.f, 2.f, 3.f});
  REQUIRE(get_each(result, [](const CtrlState &s) { return s.error; })
          == std::vector<double>{-1.f, -1.f, -1.f});
  REQUIRE(get_each(result, [](const CtrlState &s) { return s.errSum; })
          == std::vector<double>{1.f, 2.f, 3.f});
  REQUIRE(get_each(result, [](const CtrlState &s) { return s.setPt; })
          == std::vector<double>{0.f, 0.f, 0.f});
}

TEST_CASE(
    "A derivative controller with kd = 1 and setpoint zero, given a unit "
    "amplitude sawtooth should demonstrate the appropriate "
    "derivative/difference output. The error accumulation should also reflect "
    "the sawtooth signal.") {
  const auto foldable_00d = pid_algebra<double>(0.f, 0.f, 1.f);

  RealSignalVector ones{{now + 1s, 0}, {now + 2s, 1}, {now + 3s, 0}};

  CtrlState init{now, 0.f, 0.f, 0.f, 0.f};

  auto result = list_fold(foldable_00d, init, ones);

  REQUIRE(get_each(result, [](const CtrlState &s) { return s.time; })
          == std::vector<std::decay_t<decltype(now)>>{now + 1s, now + 2s,
                                                      now + 3s});
  REQUIRE(get_each(result, [](const CtrlState &s) { return s.ctrlVal; })
          == std::vector<double>{0.f, 1.f, -1.f});
  REQUIRE(get_each(result, [](const CtrlState &s) { return s.error; })
          == std::vector<double>{0.f, -1.f, 0.f});
  REQUIRE(get_each(result, [](const CtrlState &s) { return s.errSum; })
          == std::vector<double>{0.f, 1.f, 1.f});
  REQUIRE(get_each(result, [](const CtrlState &s) { return s.setPt; })
          == std::vector<double>{0.f, 0.f, 0.f});
}
