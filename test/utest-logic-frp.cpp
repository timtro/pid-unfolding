#include <sodium/sodium.h>
#include <catch/catch.hpp>
#include <iostream>

#include "../lib/pid.hpp"
#include "Plant.hpp"

using Real_t = double;

// alg : X → U → U
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

  cError.send({now + 1s, 1});
  REQUIRE(cControl.sample().value == 1.);
  REQUIRE(cControl.sample().time == now + 2s);

  cError.send({now + 2s, 1});
  REQUIRE(cControl.sample().value == 2.);
  REQUIRE(cControl.sample().time == now + 3s);
}
