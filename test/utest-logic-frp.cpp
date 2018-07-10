#include <sodium/sodium.h>
#include <catch/catch.hpp>
#include <iostream>

#include "../lib/control-frp.hpp"
#include "../lib/pid.hpp"
#include "../lib/Plant.hpp"

#include "util.hpp"

using Real_t = double;

TEST_CASE(
    "Given that plant and controller states are just timestamped doubles "
    "(i.e., signals), and an F-algebra which adds a second to the timestamps "
    "and accumulates a sum in the signal value …") {
  using CState = SignalPt<Real_t>;
  using PState = SignalPt<Real_t>;

  const auto now = chrono::steady_clock::now();

  constexpr auto inc_sum_alg = [](const PState& x, const CState& u) -> CState {
    return CState{x.time + 1s, u.value + x.value};
  };

  SECTION(
      "The sequence {(now+1s, 1), (now+2s, 3), (now+3s, 5)} pushed into the "
      "stream should produce the expected values given the algebra.") {
    const sodium::stream_sink<PState> cError;
    const CState e0 = {now, 0};

    const auto cControl = ctrl::control_scan(inc_sum_alg, cError, e0);

    auto out = std::make_shared<std::vector<CState>>();
    auto unlisten_cControl =
        cControl.listen([out](auto u) { out->push_back(u); });

    cError.send({now + 1s, 1});
    cError.send({now + 2s, 3});
    cError.send({now + 3s, 5});

    unlisten_cControl();

    REQUIRE(util::fmap([](auto x) { return x.time; }, *out)
            == std::vector{now, now + 2s, now + 3s, now + 4s});
    REQUIRE(util::fmap([](auto x) { return x.value; }, *out)
            == std::vector{0., 1., 4., 9.});
  }

  SECTION(
      "Feedback via sodium::cell::sample(). Because the result of the running "
      "sum in inc_sum_alg is fed back, the result is to double the accumulated "
      "value each time.") {
    // Why the doubling?
    //     ╭───╮  ╭───╮  ╭───╮  ╭───╮ -- feedback.
    // {1, 1 + 1, 2 + 2, 4 + 4, 8 + 8, …}
    //  ╰──╯ ╰────╯ ╰────╯ ╰────╯     -- accumulation

    const sodium::stream_sink<PState> cError;
    const CState e0 = {now, 1};
    const auto cControl = ctrl::control_scan(inc_sum_alg, cError, e0);

    auto out = std::make_shared<std::vector<CState>>();
    auto unlisten_cControl =
        cControl.listen([out](auto u) { out->push_back(u); });

    for (int k = 0; k <= 3; ++k) {
      auto ctrl = cControl.sample();
      ctrl.time -= 1s;  // Undo the increment from inc_sum_alg. Just 'cus.
      cError.send({ctrl.time, ctrl.value});
    }

    unlisten_cControl();

    REQUIRE(util::fmap([](auto x) { return x.time; }, *out)
            == std::vector{e0.time, e0.time, e0.time, e0.time, e0.time});
    REQUIRE(util::fmap([](auto x) { return x.value; }, *out)
            == std::vector{e0.value, 2., 4., 8., 16.});
  }
}
