#include <sodium/sodium.h>
#include <catch/catch.hpp>
#include <iostream>

#include "../lib/pid.hpp"
#include "Plant.hpp"

using Real_t = double;

template <typename Alg, typename U>
auto control_frp(Alg alg, const sodium::stream<U> cError, const U u0) {
  return cError.template accum<U>(u0, alg);
}

TEST_CASE(
    "Given a stream of `SignalPt<double>` for controller state and plant "
    "state, and an algebra which increments the time component and accumulates "
    "the value component, the `control_frp` should return a stream that folds "
    "the algebra over the input stream.") {
  using CState = SignalPt<Real_t>;
  using SState = SignalPt<Real_t>;

  const auto now = chrono::steady_clock::now();

  const sodium::stream_sink<SState> cError;
  const CState e0 = {now, 0};

  constexpr auto sum_alg = [](const SState& x, const CState& u) -> CState {
    return CState{x.time + 1s, u.value + x.value};
  };

  const auto cControl = control_frp(sum_alg, cError, e0);

  auto out = std::make_shared<std::vector<CState>>();
  auto unlisten_cControl =
      cControl.listen([out](auto u) { out->push_back(u); });

  cError.send({now + 1s, 1});
  cError.send({now + 2s, 3});
  cError.send({now + 3s, 5});

  unlisten_cControl();

  REQUIRE((*out)[0].value == 0.);
  REQUIRE((*out)[0].time == now);
  REQUIRE((*out)[1].value == 1.);
  REQUIRE((*out)[1].time == now + 2s);
  REQUIRE((*out)[2].value == 6.);
  REQUIRE((*out)[2].time == now + 3s);
  REQUIRE((*out)[2].value == 6.);
  REQUIRE((*out)[2].time == now + 3s);

}
