#include <cmath>
#include <iomanip>
#include <iostream>

#include <boost/numeric/odeint.hpp>
#include <catch/catch.hpp>
#include <rxcpp/rx.hpp>

#include "../include/Plant.hpp"
#include "../include/pid.hpp"
#include "../include/util/util.hpp"

#include "calculations/analytical_solutions.cpp"

#ifdef PLOT
#include "../include/plotting/gnuplot-iostream.h"
#include "../include/plotting/plot-helpers.hpp"
#endif  // PLOT

using namespace std::string_literals;

namespace ode = boost::numeric::odeint;
using CState = PIDState<>;                       // AKA U
using PState = SignalPt<std::array<double, 2>>;  // AKA X
// NB:
// PState = SignalPt<std::array<double, 2>>
//          ┌                ┐
//          │        ┌   ┐   │
//        = │  · ,   │ · │   │
//          │        │ · │   │
//          │        └   ┘   │
//          └                ┘
//            ^ Time   ^ [Position, Speed]
//
// On the other hand, sim::PState is just for Boost.odeint. It doesn't need
// time, but must be augmented with the control variable:
//
// sim::Pstate = std::array<double, 3>
//                ┌   ┐
//                │ · │  // Position
//             =  │ · │  // Speed
//                │ · │  // control variable for Boost.odeint.
//                └   ┘

constexpr double dt = 0.01;  // seconds.
constexpr auto dts = util::double_to_duration(dt);
const auto now = chrono::steady_clock::now();

constexpr double mass = 1.;
constexpr double damp = 10. / mass;
constexpr double spring = 20. / mass;
constexpr double staticForce = 1. / mass;
constexpr auto simDuration = 2s;  // seconds

// position_error : (PState, PState) → SignalPt<double>
inline SignalPt<double> position_error(const std::tuple<PState, PState>& ab) {
  auto& [a, b] = ab;
  return {std::max(a.time, b.time), a.value[0] - b.value[0]};
}

struct WorldInterface {
  const rxcpp::subjects::behavior<PState> txSubject;
  const rxcpp::observable<PState> setPoint =
      // Setpoint to x = 1, for step response.
      rxcpp::observable<>::just(PState{now, {1., 0.}});

  WorldInterface(PState x0) : txSubject(x0) {}

  void controlled_step(CState u) {
    auto x = txSubject.get_value();
    sim::PState xAugmented = {x.value[0], x.value[1], u.ctrlVal};

    if ((x.time - now) >= simDuration)
      txSubject.get_subscriber().on_completed();

    // do_step uses the second argument for both input and output.
    _stepper.do_step(_plant, xAugmented, 0, dt);
    x.time += dts;
    x.value = {xAugmented[0], xAugmented[1]};
    txSubject.get_subscriber().on_next(x);
  };

  auto get_state_observable() { return txSubject.get_observable(); }
  auto time_elapsed() { return txSubject.get_value().time - now; }

 private:
  ode::runge_kutta4<sim::PState> _stepper;
  const sim::Plant _plant = sim::Plant(staticForce, damp, spring);
};

void step_response_test(const std::string testTitle, const double Kp,
                        const double Ki, const double Kd,
                        const std::function<double(double)> expected_fn,
                        const double margin) {
  const CState u0 = {now - dts, 0., 0., 0.};
  const PState x0 = {now, {0., 0.}};
  WorldInterface worldIface(x0);

  std::vector<PState> plantStateRecord;
  worldIface.get_state_observable().subscribe(
      [&](PState x) { plantStateRecord.push_back(x); });

  // A classical confiuration for a feedback controller is illustrated as:
  //                  err   u
  //   setPoint ──➤ ⊕ ──➤ C ──➤ P ──┬──➤ plantState
  //               -↑               │
  //                │               │
  //                ╰───────────────┘
  // If you open the loop and place an interface to the imperative world, ◼,
  // at the endpoints, the controller becomes:
  //
  // setPoint  err   u
  //    ◼ ─➤ ⊕ ──➤ C ──➤ ◼
  //        -↑
  //         │ plantState
  //         ◼
  // which is in 1:1 correspondence with the code below (read backward from C)
  //
  const auto sControls =
      worldIface
          .get_state_observable()               // worldIface == ◼.
          .combine_latest(worldIface.setPoint)  // ┐ Combine state and setPt,
          .map(&position_error)                 // ┘   and then ⊕.
          .observe_on(rxcpp::identity_current_thread())
          .scan(u0, pid_algebra(Kp, Ki, Kd));  //   This is C

  sControls.subscribe(
      [&worldIface](CState u) { worldIface.controlled_step(u); });

  {
    auto simulatedPositions =
        util::fmap([](auto x) { return x.value[0]; }, plantStateRecord);
    auto theoreticalPositions = util::fmap(
        [&](auto x) { return expected_fn(util::unchrono_sec(x.time - now)); },
        plantStateRecord);

#ifdef PLOT
    const auto testData = util::fmap(
        [](const auto& x) {
          return std::make_pair(util::unchrono_sec(x.time - now), x.value[0]);
        },
        plantStateRecord);

    plot_with_tube(testTitle, testData, expected_fn, margin);
#endif  // PLOT

    REQUIRE(
        util::compareVectors(simulatedPositions, theoreticalPositions, margin));
  }
}

TEST_CASE(
    "Given system and controller parameters, simulation should reproduce "
    "analytically computed step responses to within a margin of error. See "
    "src/calculations for details. Simulations performed using FRx.") {
  SECTION("Test A (Proportional Control) FRx.") {
    constexpr double Kp = 300.;
    constexpr double Ki = 0.;
    constexpr double Kd = 0.;
    const auto title = "Test A; (Kp, Ki, Kd) = (300., 0., 0.); FRx"s;
    step_response_test(title, Kp, Ki, Kd, &analyt::test_A, 0.1);
  }

  SECTION("Test B (Proportional-Derivative Control) FRx.") {
    constexpr double Kp = 300.;
    constexpr double Ki = 0.;
    constexpr double Kd = 10.;
    const auto title = "Test B; (Kp, Ki, Kd) = (300., 0., 10.); FRx"s;
    step_response_test(title, Kp, Ki, Kd, &analyt::test_B, 0.1);
  }

  SECTION("Test C (Proportional-Integral Control) FRx.") {
    constexpr double Kp = 30.;
    constexpr double Ki = 70.;
    constexpr double Kd = 0.;
    const auto title = "Test C; (Kp, Ki, Kd) = (30., 70., 0.); FRx"s;
    step_response_test(title, Kp, Ki, Kd, &analyt::test_C, 0.1);
  }

  SECTION("Test D (Proportional-Integral-Derivative Control) FRx.") {
    constexpr double Kp = 350.;
    constexpr double Ki = 300.;
    constexpr double Kd = 50.;
    const auto title = "Test D; (Kp, Ki, Kd) = (350., 300., 50.); FRx"s;
    step_response_test(title, Kp, Ki, Kd, &analyt::test_D, 0.2);
  }
}
