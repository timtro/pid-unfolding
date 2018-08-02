#include <cmath>
#include <iomanip>
#include <iostream>

#include <boost/numeric/odeint.hpp>
#include <catch/catch.hpp>
#include <rxcpp/rx.hpp>

#include "../include/Plant.hpp"
#include "../include/control-frp.hpp"
#include "../include/pid.hpp"
#include "../include/util/util-frp.hpp"
#include "../include/util/util.hpp"

#include "calculations/analytical_solutions.cpp"

#ifdef PLOT
#include "../include/plotting/gnuplot-iostream.h"
#include "../include/plotting/plot-helpers.hpp"
#endif  // PLOT

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

inline SignalPt<double> position_error(const std::tuple<PState, PState>& ab) {
  auto& [a, b] = ab;
  return {std::max(a.time, b.time), a.value[0] - b.value[0]};
}

struct WorldInterface {
  const rxcpp::subjects::behavior<PState> txSubject;
  const rxcpp::observable<PState> setPoint =
      // Setpoint to x = 1, for step response.
      rxcpp::observable<>::just(PState{now, {1., 0.}});

  WorldInterface(PState x0) : txSubject(x0), _x(x0) {}

  void controlled_step(CState u) {
    sim::PState xAugmented = {_x.value[0], _x.value[1], u.ctrlVal};

    // do_step uses the second argument for both input and output.
    _stepper.do_step(_plant, xAugmented, 0, dt);
    _x.time += dts;
    _x.value = {xAugmented[0], xAugmented[1]};
    txSubject.get_subscriber().on_next(_x);
  };

  auto get_state_observable() { return txSubject.get_observable(); }
  auto time_elapsed() { return _x.time - now; }

 private:
  PState _x;
  ode::runge_kutta4<sim::PState> _stepper;
  const sim::Plant _plant = sim::Plant(staticForce, damp, spring);
};

TEST_CASE(
    "Given system and controller parameters, simulation should reproduce "
    "analytically computed step responses to within a margin of error. See "
    "src/calculations for details. Simulations performed using FRP.") {
  const CState u0 = {now - dts, 0., 0., 0.};
  const PState x0 = {now, {0., 0.}};
  WorldInterface worldIface(x0);

  SECTION("Test A (Proportional Control) FRP.") {
    constexpr double Kp = 300.;
    constexpr double Ki = 0.;
    constexpr double Kd = 0.;

    std::vector<PState> plantStateRecord;
    worldIface.get_state_observable().subscribe(
        [&](PState x) { plantStateRecord.push_back(x); });

    CState controllerState;

    // A classical confiuration for a feedback controller is illustrated as:
    //
    //   setPoint ──➤ ⊕ ──➤ C ──➤ P ──┬──➤ plantState
    //               -↑               │
    //                │               │
    //                ╰───────────────┘
    // If you open the loop and place an interface to the imperative world, ◼,
    // at the endpoints, the controller becomes:
    //
    //    ◼ ─➤ ⊕ ──➤ C ──➤ ◼
    //        -↑
    //         │
    //         ◼
    // which is in 1:1 correspondence with the code below (read backward from C)
    //
    const auto sControls =
        worldIface.get_state_observable()          // worldIface == ◼.
            .combine_latest(worldIface.setPoint)   // ┐ Combine state and setPt,
            .map(&position_error)                  // ┘   and then ⊕.
            .scan(u0, pid_algebra(Kp, Ki, Kd));//   This is C

    sControls.subscribe([&controllerState](CState u) { controllerState = u; });

    // Imperative loop where the FRx meets the real world.
    while (worldIface.time_elapsed() < simDuration) {
      worldIface.controlled_step(controllerState);
    }

    {
      constexpr double margin = 0.1;
      auto simulatedPositions =
          util::fmap([](auto x) { return x.value[0]; }, plantStateRecord);
      auto theoreticalPositions = util::fmap(
          [](auto x) {
            return analyt::test_A(util::unchrono_sec(x.time - now));
          },
          plantStateRecord);

#ifdef PLOT
      const auto testData = util::fmap(
          [](const auto& x) {
            return std::make_pair(util::unchrono_sec(x.time - now), x.value[0]);
          },
          plantStateRecord);

      plot_with_tube("Test A, (Kp, Ki, Kd) = (300, 0, 0).", testData,
                     &analyt::test_A, margin);
#endif  // PLOT

      REQUIRE(util::compareVectors(simulatedPositions, theoreticalPositions,
                                   margin));
    }
  }

  //   SECTION("Test B (Proportional-Derivative Control) anamorphism.") {
  //     constexpr double Kp = 300.;
  //     constexpr double Ki = 0.;
  //     constexpr double Kd = 10.;
  //
  //     const sodium::stream_sink<SignalPt<PState>> sPlantState;
  //
  //     auto cControl = ctrl::control_frp(
  //         pid_algebra(Kp, Ki, Kd), &to_error, setPoint,
  //         static_cast<sodium::stream<SignalPt<PState>>>(sPlantState), u0);
  //
  //     auto [sPlantState_unlisten, plantRecord] =
  //     util::make_listener(sPlantState);
  //     {
  //       PState x = {0., 0., 0.};
  //       sPlantState.send({now, x});
  //       for (int k = 1; k < simDuration / dt; ++k) {
  //         x[2] = cControl.sample().ctrlVal;
  //         stepper.do_step(plant, x, 0, dt);
  //         sPlantState.send({now + k * dts, x});
  //       }
  //     }
  //     sPlantState_unlisten();
  //
  //     {
  //       constexpr double margin = 0.1;
  //       auto simulatedPositions =
  //           util::fmap([](auto x) { return x.value[0]; }, *plantRecord);
  //       auto theoreticalPositions = util::fmap(
  //           [](auto x) {
  //             return analyt::test_B(util::unchrono_sec(x.time - now));
  //           },
  //           *plantRecord);
  //
  // #ifdef PLOT
  //       const auto testData = util::fmap(
  //           [](const auto& x) {
  //             return std::make_pair(util::unchrono_sec(x.time - now),
  //             x.value[0]);
  //           },
  //           *plantRecord);
  //
  //       plot_with_tube("Test B, (Kp, Ki, Kd) = (300, 0, 10).", testData,
  //                      &analyt::test_B, margin);
  // #endif  // PLOT
  //
  //       REQUIRE(util::compareVectors(simulatedPositions,
  //       theoreticalPositions,
  //                                    margin));
  //     }
  //   }
  //
  //   SECTION("Test C (Proportional-Integral Control) anamorphism.") {
  //     constexpr double Kp = 30.;
  //     constexpr double Ki = 70.;
  //     constexpr double Kd = 0.;
  //
  //     // Setpoint to x = 1, for step response.
  //     const sodium::cell<PState> setPoint({1., 0., 0.});
  //     const sodium::stream_sink<SignalPt<PState>> sPlantState;
  //
  //     const CState u0 = {now - dts, 0., 0., 0.};
  //     auto cControl = ctrl::control_frp(
  //         pid_algebra(Kp, Ki, Kd), &to_error, setPoint,
  //         static_cast<sodium::stream<SignalPt<PState>>>(sPlantState), u0);
  //
  //     auto [sPlantState_unlisten, plantRecord] =
  //     util::make_listener(sPlantState);
  //     {
  //       PState x = {0., 0., 0.};
  //       sPlantState.send({now, x});
  //       for (int k = 1; k < simDuration / dt; ++k) {
  //         x[2] = cControl.sample().ctrlVal;
  //         stepper.do_step(plant, x, 0, dt);
  //         sPlantState.send({now + k * dts, x});
  //       }
  //     }
  //     sPlantState_unlisten();
  //
  //     {
  //       constexpr double margin = 0.1;
  //
  //       auto simulatedPositions =
  //           util::fmap([](auto x) { return x.value[0]; }, *plantRecord);
  //       auto theoreticalPositions = util::fmap(
  //           [](auto x) {
  //             return analyt::test_C(util::unchrono_sec(x.time - now));
  //           },
  //           *plantRecord);
  //
  // #ifdef PLOT
  //       const auto testData = util::fmap(
  //           [](const auto& x) {
  //             return std::make_pair(util::unchrono_sec(x.time - now),
  //             x.value[0]);
  //           },
  //           *plantRecord);
  //
  //       plot_with_tube("Test C, (Kp, Ki, Kd) = (30, 70, 0).", testData,
  //                      &analyt::test_C, margin);
  // #endif  // PLOT
  //
  //       REQUIRE(util::compareVectors(simulatedPositions,
  //       theoreticalPositions,
  //                                    margin));
  //     }
  //   }
  //
  //   SECTION("Test D (Proportional-Integral-Derivative Control) anamorphism.")
  //   {
  //     constexpr double Kp = 350.;
  //     constexpr double Ki = 300.;
  //     constexpr double Kd = 50.;
  //
  //     // Setpoint to x = 1, for step response.
  //     const sodium::cell<PState> setPoint({1., 0., 0.});
  //     const sodium::stream_sink<SignalPt<PState>> sPlantState;
  //
  //     const CState u0 = {now - dts, 0., 0., 0.};
  //     auto cControl = ctrl::control_frp(
  //         pid_algebra(Kp, Ki, Kd), &to_error, setPoint,
  //         static_cast<sodium::stream<SignalPt<PState>>>(sPlantState), u0);
  //
  //     auto [sPlantState_unlisten, plantRecord] =
  //     util::make_listener(sPlantState);
  //     {
  //       PState x = {0., 0., 0.};
  //       sPlantState.send({now, x});
  //       for (int k = 1; k < simDuration / dt; ++k) {
  //         x[2] = cControl.sample().ctrlVal;
  //         stepper.do_step(plant, x, 0, dt);
  //         sPlantState.send({now + k * dts, x});
  //       }
  //     }
  //     sPlantState_unlisten();
  //
  //     {
  //       constexpr double margin = 0.2;
  //
  //       auto simulatedPositions =
  //           util::fmap([](auto x) { return x.value[0]; }, *plantRecord);
  //       auto theoreticalPositions = util::fmap(
  //           [](auto x) {
  //             return analyt::test_D(util::unchrono_sec(x.time - now));
  //           },
  //           *plantRecord);
  //
  // #ifdef PLOT
  //       const auto testData = util::fmap(
  //           [](const auto& x) {
  //             return std::make_pair(util::unchrono_sec(x.time - now),
  //             x.value[0]);
  //           },
  //           *plantRecord);
  //
  //       plot_with_tube("Test B, (Kp, Ki, Kd) = (350, 300, 50).", testData,
  //                      &analyt::test_D, margin);
  // #endif  // PLOT
  //
  //       REQUIRE(util::compareVectors(simulatedPositions,
  //       theoreticalPositions,
  //                                    margin));
  //     }
  //   }
}
