#include <boost/numeric/odeint.hpp>
#include <catch/catch.hpp>
#include <range/v3/all.hpp>

#include "../include/Plant.hpp"
#include "../include/control-frp.hpp"
#include "../include/pid.hpp"
#include "../include/util/util.hpp"

#include "calculations/analytical_solutions.cpp"

#ifdef PLOT
#include "../include/plotting/gnuplot-iostream.h"
#include "../include/plotting/plot-helpers.hpp"
#endif  // PLOT

using namespace ranges;
namespace ode = boost::numeric::odeint;
using CState = PIDState<double>;
using PState = sim::PState;

constexpr double dt = 0.01;  // seconds.
constexpr auto dts = util::double_to_duration(dt);

const auto now = chrono::steady_clock::now();

constexpr double mass = 1.;
constexpr double damp = 10. / mass;
constexpr double spring = 20. / mass;
constexpr double staticForce = 1. / mass;
constexpr double simTime = 2;  // seconds
const sim::Plant plant(staticForce, damp, spring);

ode::runge_kutta4<PState> stepper;

inline SignalPt<double> to_error(const SignalPt<PState>& a, const PState& b) {
  return {a.time, a.value[0] - b[0]};
}

const std::vector<chrono::time_point<chrono::steady_clock>> ts =
    view::iota(1) | view::transform([](auto k) { return now + k * dts; })
    | view::take_while([](auto t) { return t > (now + 2s); });

// TEST_CASE(
//     "Test A (Proportional Control)—Reproduct known (analytical) result, "
//     "borrowed from\n"
//     "http://ctms.engin.umich.edu/CTMS/"
//     "index.php?example=Introduction&section=ControlPID\n"
//     "A damped-driven harmonic oscillator under the influence of a "
//     "P-controller, both with parameters defined in the test, should produce
//     a" "step-response within a margin of the analytical solution.",
//     "[Test 1], [P-controller]") {
//   constexpr double Kp = 300.;
//   constexpr double Ki = 0.;
//   constexpr double Kd = 0.;
//
//   const auto pid_controller = pid_algebra(Kp, Ki, Kd);
//
//   const CState u0 = {now - dts, 0., 0., 0.};
//   const PState x0 = {0., 0., 0.};
//
//   const auto stepResponseCoalgebra = [&pid_controller](auto xu)
//       -> std::optional<std::pair<std::pair<SignalPt<PState>, CState>,
//                                  std::pair<SignalPt<PState>, CState>>> {
//     const auto [tx, u] = xu;
//
//     const auto thisTime = tx.time;
//     const auto x = tx.value;
//     assert(thisTime == u.time);
//     const auto nextTime = thisTime + dts;
//
//     CState nextControllerState = pid_controller();
//     PState nextPlantState = sim::do_step_with(plant, stepper, x);
//
//     return {{}, {}};
//   };
//
//   {
//     constexpr double margin = 0.1;
//     auto frpPositions =
//         util::fmap([](auto x) { return x.value[0]; }, *plantRecord);
//     auto realPositions = util::fmap(
//         [](auto x) { return analyt::test_A(util::unchrono_sec(x.time - now));
//         }, *plantRecord);
//
// #ifdef PLOT
//
//     const auto testData = util::fmap(
//         [](const auto& x) {
//           return std::make_pair(util::unchrono_sec(x.time - now),
//           x.value[0]);
//         },
//         *plantRecord);
//
//     plot_with_tube("Test A, (Kp, Ki, Kd) = (300, 0, 0).", testData,
//                    &analyt::test_A, margin);
//
// #endif  // PLOT
//
//     REQUIRE(util::compareVectors(frpPositions, realPositions, margin));
//   }
// }

// TEST_CASE(
//     "Test B (Proportional-Derivative Control)—Reproduct known (analytical) "
//     "result, borrowed from\n"
//     "http://ctms.engin.umich.edu/CTMS/"
//     "index.php?example=Introduction&section=ControlPID\n"
//     "A damped-driven harmonic oscillator under the influence of a "
//     "PD-controller, both with parameters defined in the test, should produce
//     a " "step-response within a margin of the analytical solution.",
//     "[Test 1], [P-controller]") {
//   constexpr double Kp = 300.;
//   constexpr double Ki = 0.;
//   constexpr double Kd = 10.;
//
//   // Setpoint to x = 1, for step response.
//   const sodium::cell<PState> setPoint({1., 0., 0.});
//   const sodium::stream_sink<SignalPt<PState>> sPlantState;
//
//   const CState u0 = {now - dts, 0., 0., 0.};
//   auto cControl = ctrl::control_frp(
//       pid_algebra(Kp, Ki, Kd), &to_error, setPoint,
//       static_cast<sodium::stream<SignalPt<PState>>>(sPlantState), u0);
//
//   auto [sPlantState_unlisten, plantRecord] =
//   util::make_listener(sPlantState);
//   {
//     PState x = {0., 0., 0.};
//     sPlantState.send({now, x});
//     for (int k = 1; k < simTime / dt; ++k) {
//       x[2] = cControl.sample().ctrlVal;
//       stepper.do_step(plant, x, 0, dt);
//       sPlantState.send({now + k * dts, x});
//     }
//   }
//   sPlantState_unlisten();
//
//   {
//     constexpr double margin = 0.1;
//     auto frpPositions =
//         util::fmap([](auto x) { return x.value[0]; }, *plantRecord);
//     auto realPositions = util::fmap(
//         [](auto x) { return analyt::test_B(util::unchrono_sec(x.time - now));
//         }, *plantRecord);
//
// #ifdef PLOT
//     const auto testData = util::fmap(
//         [](const auto& x) {
//           return std::make_pair(util::unchrono_sec(x.time - now),
//           x.value[0]);
//         },
//         *plantRecord);
//
//     plot_with_tube("Test B, (Kp, Ki, Kd) = (300, 0, 10).", testData,
//                    &analyt::test_B, margin);
// #endif  // PLOT
//
//     REQUIRE(util::compareVectors(frpPositions, realPositions, margin));
//   }
// }
//
// TEST_CASE(
//     "Test C (Proportional-Integral Control)—Reproduct known (analytical) "
//     "result, borrowed from\n"
//     "http://ctms.engin.umich.edu/CTMS/"
//     "index.php?example=Introduction&section=ControlPID\n"
//     "A damped-driven harmonic oscillator under the influence of a "
//     "PI-controller, both with parameters defined in the test, should produce
//     a " "step-response within a margin of the analytical solution.",
//     "[Test 1], [P-controller]") {
//   constexpr double Kp = 30.;
//   constexpr double Ki = 70.;
//   constexpr double Kd = 0.;
//
//   // Setpoint to x = 1, for step response.
//   const sodium::cell<PState> setPoint({1., 0., 0.});
//   const sodium::stream_sink<SignalPt<PState>> sPlantState;
//
//   const CState u0 = {now - dts, 0., 0., 0.};
//   auto cControl = ctrl::control_frp(
//       pid_algebra(Kp, Ki, Kd), &to_error, setPoint,
//       static_cast<sodium::stream<SignalPt<PState>>>(sPlantState), u0);
//
//   auto [sPlantState_unlisten, plantRecord] =
//   util::make_listener(sPlantState);
//   {
//     PState x = {0., 0., 0.};
//     sPlantState.send({now, x});
//     for (int k = 1; k < simTime / dt; ++k) {
//       x[2] = cControl.sample().ctrlVal;
//       stepper.do_step(plant, x, 0, dt);
//       sPlantState.send({now + k * dts, x});
//     }
//   }
//   sPlantState_unlisten();
//
//   {
//     constexpr double margin = 0.1;
//
//     auto frpPositions =
//         util::fmap([](auto x) { return x.value[0]; }, *plantRecord);
//     auto realPositions = util::fmap(
//         [](auto x) { return analyt::test_C(util::unchrono_sec(x.time - now));
//         }, *plantRecord);
//
// #ifdef PLOT
//     const auto testData = util::fmap(
//         [](const auto& x) {
//           return std::make_pair(util::unchrono_sec(x.time - now),
//           x.value[0]);
//         },
//         *plantRecord);
//
//     plot_with_tube("Test C, (Kp, Ki, Kd) = (30, 70, 0).", testData,
//                    &analyt::test_C, margin);
// #endif  // PLOT
//
//     REQUIRE(util::compareVectors(frpPositions, realPositions, margin));
//   }
// }
//
// TEST_CASE(
//     "Test D (Proportional-Integral-Derivative Control)—Reproduct known "
//     "(analytical) result, borrowed from\n"
//     "http://ctms.engin.umich.edu/CTMS/"
//     "index.php?example=Introduction&section=ControlPID\n"
//     "A damped-driven harmonic oscillator under the influence of a "
//     "PID-controller, both with parameters defined in the test, should produce
//     " "a step-response within a margin of the analytical solution.",
//     "[Test D], [P-controller]") {
//   constexpr double Kp = 350.;
//   constexpr double Ki = 300.;
//   constexpr double Kd = 50.;
//
//   // Setpoint to x = 1, for step response.
//   const sodium::cell<PState> setPoint({1., 0., 0.});
//   const sodium::stream_sink<SignalPt<PState>> sPlantState;
//
//   const CState u0 = {now - dts, 0., 0., 0.};
//   auto cControl = ctrl::control_frp(
//       pid_algebra(Kp, Ki, Kd), &to_error, setPoint,
//       static_cast<sodium::stream<SignalPt<PState>>>(sPlantState), u0);
//
//   auto [sPlantState_unlisten, plantRecord] =
//   util::make_listener(sPlantState);
//   {
//     PState x = {0., 0., 0.};
//     sPlantState.send({now, x});
//     for (int k = 1; k < simTime / dt; ++k) {
//       x[2] = cControl.sample().ctrlVal;
//       stepper.do_step(plant, x, 0, dt);
//       sPlantState.send({now + k * dts, x});
//     }
//   }
//   sPlantState_unlisten();
//
//   {
//     constexpr double margin = 0.2;
//
//     auto frpPositions =
//         util::fmap([](auto x) { return x.value[0]; }, *plantRecord);
//     auto realPositions = util::fmap(
//         [](auto x) { return analyt::test_D(util::unchrono_sec(x.time - now));
//         }, *plantRecord);
//
// #ifdef PLOT
//     const auto testData = util::fmap(
//         [](const auto& x) {
//           return std::make_pair(util::unchrono_sec(x.time - now),
//           x.value[0]);
//         },
//         *plantRecord);
//
//     plot_with_tube("Test B, (Kp, Ki, Kd) = (350, 300, 50).", testData,
//                    &analyt::test_D, margin);
// #endif  // PLOT
//
//     REQUIRE(util::compareVectors(frpPositions, realPositions, margin));
//   }
// }
