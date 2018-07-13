#include <boost/numeric/odeint.hpp>
#include <catch/catch.hpp>
#include <range/v3/all.hpp>

#include "../include/Plant.hpp"
#include "../include/control-frp.hpp"
#include "../include/pid.hpp"
#include "../include/util/util-sim.hpp"
#include "../include/util/util.hpp"

#include "calculations/analytical_solutions.cpp"

#ifdef PLOT
#include "../include/plotting/gnuplot-iostream.h"
#include "../include/plotting/plot-helpers.hpp"
#endif  // PLOT

using namespace ranges;
namespace ode = boost::numeric::odeint;

using CState = PIDState<>;                       // AKA U
using PState = SignalPt<std::array<double, 2>>;  // AKA X
// NB:
// PState = SignalPt<std::array<double, 2>>
//          ┌                ┐
//          │        ┌   ┐   │
//        = │  · ,   │ · │   │
//          │        │ · │   │
//          └        └   ┘   ┘
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
constexpr double simTime = 2;  // seconds
const sim::Plant plant(staticForce, damp, spring);
ode::runge_kutta4<sim::PState> stepper;

const CState u0 = {now - dts, 0., 0., 0.};
const PState x0 = {now, {0., 0.}};

// controlled_step : (X, U) → X  = (PState, CState) → PState
const auto controlled_step = [](const PState& x, const CState& u) -> PState {
  // Remember: sim::PState is the timeless, augmented PState, not X. (See NB.)
  sim::PState xOld = {x.value[0], x.value[1], u.ctrlVal};
  sim::PState xNew = util::do_step_with(plant, stepper, dt, xOld);

  return {x.time + dts, {xNew[0], xNew[1]}};
};

// step_response_coalg : (X, U) → optional<((X, U), (X, U))>
const auto make_step_response_coalg = [](auto controller) {
  return
      [controller](const std::pair<PState, CState>& xu)
          -> std::optional<
              std::pair<std::pair<PState, CState>, std::pair<PState, CState>>> {
        const auto [x, u] = xu;
        constexpr double positionSetpoint = 1.;

        if (x.time > now + 2s) return {};

        const SignalPt<double> error = {x.time, x.value[0] - positionSetpoint};
        CState uNew = controller(error, u);
        PState xNew = controlled_step(x, uNew);

        return {{{xNew, uNew}, {x, u}}};
      };
};

TEST_CASE(
    "Given system and controller parameters, simulation should reproduce "
    "analytically computed step responses to within a margin of error. See "
    "src/calculations for details. Simulations performed using anamorphism.") {
  SECTION("Test A (Proportional Control) anamorphism.") {
    constexpr double Kp = 300.;
    constexpr double Ki = 0.;
    constexpr double Kd = 0.;

    auto step_response_coalg =
        make_step_response_coalg(pid_algebra(Kp, Ki, Kd));

    // result : std::vector<std::pair<PState, CState>>
    auto result = util::unfold(step_response_coalg, std::pair{x0, u0});

    {
      constexpr double margin = 0.1;
      auto simulatedPositions =
          util::fmap([](auto xu) { return xu.first.value[0]; }, result);
      auto theoreticalPositions = util::fmap(
          [](auto xu) {
            return analyt::test_A(util::unchrono_sec(xu.first.time - now));
          },
          result);

#ifdef PLOT

      const auto testData = util::fmap(
          [](const auto& xu) {
            return std::make_pair(util::unchrono_sec(xu.first.time - now),
                                  xu.first.value[0]);
          },
          result);

      plot_with_tube("Test A, (Kp, Ki, Kd) = (300, 0, 0).", testData,
                     &analyt::test_A, margin);

#endif  // PLOT

      REQUIRE(util::compareVectors(simulatedPositions, theoreticalPositions,
                                   margin));
    }
  }

  SECTION("Test B (Proportional-Derivative Control) anamorphism.") {
    constexpr double Kp = 300.;
    constexpr double Ki = 0.;
    constexpr double Kd = 10.;

    auto step_response_coalg =
        make_step_response_coalg(pid_algebra(Kp, Ki, Kd));

    // result : std::vector<std::pair<PState, CState>>
    auto result = util::unfold(step_response_coalg, std::pair{x0, u0});

    {
      constexpr double margin = 0.1;
      auto simulatedPositions =
          util::fmap([](auto xu) { return xu.first.value[0]; }, result);
      auto theoreticalPositions = util::fmap(
          [](auto xu) {
            return analyt::test_B(util::unchrono_sec(xu.first.time - now));
          },
          result);

#ifdef PLOT

      const auto testData = util::fmap(
          [](const auto& xu) {
            return std::make_pair(util::unchrono_sec(xu.first.time - now),
                                  xu.first.value[0]);
          },
          result);

      plot_with_tube("Test B, (Kp, Ki, Kd) = (300, 0, 10).", testData,
                     &analyt::test_B, margin);

#endif  // PLOT

      REQUIRE(util::compareVectors(simulatedPositions, theoreticalPositions,
                                   margin));
    }
  }

  SECTION("Test C (Proportional-Integral Control) anamorphism.") {
    constexpr double Kp = 30.;
    constexpr double Ki = 70.;
    constexpr double Kd = 0.;

    auto step_response_coalg =
        make_step_response_coalg(pid_algebra(Kp, Ki, Kd));

    // result : std::vector<std::pair<PState, CState>>
    auto result = util::unfold(step_response_coalg, std::pair{x0, u0});

    {
      constexpr double margin = 0.1;
      auto simulatedPositions =
          util::fmap([](auto xu) { return xu.first.value[0]; }, result);
      auto theoreticalPositions = util::fmap(
          [](auto xu) {
            return analyt::test_C(util::unchrono_sec(xu.first.time - now));
          },
          result);

#ifdef PLOT

      const auto testData = util::fmap(
          [](const auto& xu) {
            return std::make_pair(util::unchrono_sec(xu.first.time - now),
                                  xu.first.value[0]);
          },
          result);

      plot_with_tube("Test C, (Kp, Ki, Kd) = (30, 70, 0).", testData,
                     &analyt::test_C, margin);

#endif  // PLOT

      REQUIRE(util::compareVectors(simulatedPositions, theoreticalPositions,
                                   margin));
    }
  }

  SECTION("Test D (Proportional-Integral-Derivative Control) anamorphism.") {
    constexpr double Kp = 350.;
    constexpr double Ki = 300.;
    constexpr double Kd = 50.;

    auto step_response_coalg =
        make_step_response_coalg(pid_algebra(Kp, Ki, Kd));

    // result : std::vector<std::pair<PState, CState>>
    auto result = util::unfold(step_response_coalg, std::pair{x0, u0});

    {
      constexpr double margin = 0.2;
      auto simulatedPositions =
          util::fmap([](auto xu) { return xu.first.value[0]; }, result);
      auto theoreticalPositions = util::fmap(
          [](auto xu) {
            return analyt::test_D(util::unchrono_sec(xu.first.time - now));
          },
          result);

#ifdef PLOT

      const auto testData = util::fmap(
          [](const auto& xu) {
            return std::make_pair(util::unchrono_sec(xu.first.time - now),
                                  xu.first.value[0]);
          },
          result);

      plot_with_tube("Test D, (Kp, Ki, Kd) = (350, 300, 50).", testData,
                     &analyt::test_D, margin);

#endif  // PLOT

      REQUIRE(util::compareVectors(simulatedPositions, theoreticalPositions,
                                   margin));
    }
  }
}
