#include <cmath>
#include <iomanip>
#include <iostream>

#include <sodium/sodium.h>
#include <boost/format.hpp>
#include <boost/hana/functional/curry.hpp>
#include <boost/numeric/odeint.hpp>
#include <catch/catch.hpp>

#include "../lib/Plant.hpp"
#include "../lib/control-frp.hpp"
#include "../lib/pid.hpp"

#include "test-util.hpp"

#ifdef PLOT
#include "gnuplot-iostream.h"
#endif  // PLOT

namespace ode = boost::numeric::odeint;
using boost::hana::curry;
using Real_t = double;
using CState = PIDState<Real_t>;

constexpr double dt = 0.01;  // seconds.
constexpr auto dts = util::double_to_duration(dt);

inline SignalPt<double> to_error(const SignalPt<PState>& a, const PState& b) {
  return {a.time, a.value[0] - b[0]};
}

const auto now = chrono::steady_clock::now();

constexpr Real_t mass = 1.;
constexpr Real_t damp = 10. / mass;
constexpr Real_t spring = 20. / mass;
constexpr Real_t staticForce = -1. / mass;
constexpr Real_t simTime = 2;  // seconds
const sim::Plant plant(staticForce, damp, spring);

TEST_CASE(
    "Test A (Proportional Control)—Reproduct known (analytical) result, "
    "borrowed from\n"
    "http://ctms.engin.umich.edu/CTMS/"
    "index.php?example=Introduction&section=ControlPID\n"
    "A damped-driven harmonic oscillator under the influence of a "
    "P-controller, both with parameters defined in the test, should produce a "
    "step-response within a margin of the analytical solution.",
    "[Test 1], [P-controller]") {
  constexpr Real_t Kp = 300.;
  constexpr Real_t Ki = 0.;
  constexpr Real_t Kd = 0.;

  // Setpoint to x = 1, for step response.
  const sodium::cell<PState> setPoint({1., 0., 0.});
  const sodium::stream_sink<SignalPt<PState>> sPlantState;
  const auto sError = sPlantState.snapshot(setPoint, &to_error);

  const CState u0 = {now - dts, 0., 0., 0.};
  // const auto cControl = ctrl::control_frp(pid_algebra(Kp, Ki, Kd), sError,
  // u0); Or equivalently,
  const auto cControl = sError.accum<CState>(u0, pid_algebra(Kp, Ki, Kd));

  ode::runge_kutta4<PState> stepper;
  // do_step : Plant → Stepper → double → PState → PState
  //                             ^ dt     ^ x[k]   ^ x[k+1]
  constexpr auto do_step =
      curry<4>(sim::do_step_with<ode::runge_kutta4<PState>>);

  // plant_step : PState → PState
  const auto plant_step = do_step(plant, stepper, dt);

  // ---

  auto [sPlantState_unlisten, plantRecord] = util::make_listener(sPlantState);

  {
    PState x = {0., 0., 0.};
    sPlantState.send({now, x});
    for (int k = 1; k < simTime / dt; ++k) {
      x[2] = cControl.sample().ctrlVal;
      stepper.do_step(plant, x, 0, dt);
      sPlantState.send({now + k * dts, x});
    }
  }

  sPlantState_unlisten();

  {
    constexpr double margin = 0.1;
    constexpr auto sys = [](double t) {
      return std::exp(-5. * t)
                 * (-(150. * std::sin(std::sqrt(285.) * t))
                        / (31. * std::sqrt(285.))
                    - (30. * std::cos(std::sqrt(285.) * t)) / 31.)
             + 30. / 31.;
    };

    auto frpPositions =
        util::fmap([](auto x) { return x.value[0]; }, *plantRecord);
    auto realPositions = util::fmap(
        [&sys](auto x) { return sys(util::unchrono_sec(x.time - now)); },
        *plantRecord);

#ifdef PLOT
    Gnuplot gp;
    gp << "set title 'Test A: pid(" << boost::format("%.2f") % Kp << ", "
       << boost::format("%.2f") % Ki << ", " << boost::format("%.2f") % Kd
       << ")\n";
    gp << "plot '-' u 1:2:3 title 'Acceptable margin: analytical ±"
       << boost::format("%.3f") % margin
       << "' w filledcu fs solid fc rgb '#6699FF55', '-' u 1:2 title 'Test "
          "result' w l\n";
    auto range = util::fmap(
        [&sys](auto x) {
          const double t = util::unchrono_sec(x.time - now);
          return std::make_tuple(t, sys(t) + margin, sys(t) - margin);
        },
        *plantRecord);

    gp.send1d(range);
    gp.send1d(util::fmap(
        [](const auto& x) {
          return std::make_pair(util::unchrono_sec(x.time - now), x.value[0]);
        },
        *plantRecord));
    gp << "pause mouse key\nwhile (MOUSE_CHAR ne 'q') { pause mouse "
          "key; }\n";
#endif  // PLOT

    REQUIRE(util::compareVectors(frpPositions, realPositions, margin));
  }
}

TEST_CASE(
    "Test B (Proportional-Derivative Control)—Reproduct known (analytical) "
    "result, borrowed from\n"
    "http://ctms.engin.umich.edu/CTMS/"
    "index.php?example=Introduction&section=ControlPID\n"
    "A damped-driven harmonic oscillator under the influence of a "
    "PD-controller, both with parameters defined in the test, should produce a "
    "step-response within a margin of the analytical solution.",
    "[Test 1], [P-controller]") {
  constexpr Real_t Kp = 300.;
  constexpr Real_t Ki = 0.;
  constexpr Real_t Kd = 10.;

  // Setpoint to x = 1, for step response.
  const sodium::cell<PState> setPoint({1., 0., 0.});
  const sodium::stream_sink<SignalPt<PState>> sPlantState;
  const auto sError = sPlantState.snapshot(setPoint, &to_error);

  const CState u0 = {now - dts, 0., 0., 0.};
  // const auto cControl = ctrl::control_frp(pid_algebra(Kp, Ki, Kd), sError,
  // u0); Or equivalently,
  const auto cControl = sError.accum<CState>(u0, pid_algebra(Kp, Ki, Kd));

  ode::runge_kutta4<PState> stepper;
  // do_step : Plant → Stepper → double → PState → PState
  //                             ^ dt     ^ x[k]   ^ x[k+1]
  constexpr auto do_step =
      curry<4>(sim::do_step_with<ode::runge_kutta4<PState>>);

  // plant_step : PState → PState
  const auto plant_step = do_step(plant, stepper, dt);

  // ---

  auto [sPlantState_unlisten, plantRecord] = util::make_listener(sPlantState);

  {
    PState x = {0., 0., 0.};
    sPlantState.send({now, x});
    for (int k = 1; k < simTime / dt; ++k) {
      x[2] = cControl.sample().ctrlVal;
      stepper.do_step(plant, x, 0, dt);
      sPlantState.send({now + k * dts, x});
    }
  }

  sPlantState_unlisten();

  {
    constexpr double margin = 0.1;
    constexpr auto sys = [](double t) {
      return 0.10118 * std::exp(-11.145 * t) * std::sin(15.511 * t)
             - 1.0944 * std::exp(-11.145 * t) * std::cos(15.511 * t)
             + 0.15697 * std::exp(-87.708 * t) + 0.93749;
    };

    auto frpPositions =
        util::fmap([](auto x) { return x.value[0]; }, *plantRecord);
    auto realPositions = util::fmap(
        [&sys](auto x) { return sys(util::unchrono_sec(x.time - now)); },
        *plantRecord);

#ifdef PLOT
    Gnuplot gp;
    gp << "set title 'Test B: pid(" << boost::format("%.2f") % Kp << ", "
       << boost::format("%.2f") % Ki << ", " << boost::format("%.2f") % Kd
       << ")\n";
    gp << "plot '-' u 1:2:3 title 'Acceptable margin: analytical ±"
       << boost::format("%.3f") % margin
       << "' w filledcu fs solid fc rgb '#6699FF55', '-' u 1:2 title 'Test "
          "result' w l\n";
    auto range = util::fmap(
        [&sys](auto x) {
          const double t = util::unchrono_sec(x.time - now);
          return std::make_tuple(t, sys(t) + margin, sys(t) - margin);
        },
        *plantRecord);

    gp.send1d(range);
    gp.send1d(util::fmap(
        [](const auto& x) {
          return std::make_pair(util::unchrono_sec(x.time - now), x.value[0]);
        },
        *plantRecord));
    gp << "pause mouse key\nwhile (MOUSE_CHAR ne 'q') { pause mouse "
          "key; }\n";
#endif  // PLOT

    REQUIRE(util::compareVectors(frpPositions, realPositions, margin));
  }
}

TEST_CASE(
    "Test C (Proportional-Integral Control)—Reproduct known (analytical) "
    "result, borrowed from\n"
    "http://ctms.engin.umich.edu/CTMS/"
    "index.php?example=Introduction&section=ControlPID\n"
    "A damped-driven harmonic oscillator under the influence of a "
    "PI-controller, both with parameters defined in the test, should produce a "
    "step-response within a margin of the analytical solution.",
    "[Test 1], [P-controller]") {
  constexpr Real_t Kp = 30.;
  constexpr Real_t Ki = 70.;
  constexpr Real_t Kd = 0.;

  // Setpoint to x = 1, for step response.
  const sodium::cell<PState> setPoint({1., 0., 0.});
  const sodium::stream_sink<SignalPt<PState>> sPlantState;
  const auto sError = sPlantState.snapshot(setPoint, &to_error);

  const CState u0 = {now - dts, 0., 0., 0.};
  // const auto cControl = ctrl::control_frp(pid_algebra(Kp, Ki, Kd), sError,
  // u0); Or equivalently,
  const auto cControl = sError.accum<CState>(u0, pid_algebra(Kp, Ki, Kd));

  ode::runge_kutta4<PState> stepper;
  // do_step : Plant → Stepper → double → PState → PState
  //                             ^ dt     ^ x[k]   ^ x[k+1]
  constexpr auto do_step =
      curry<4>(sim::do_step_with<ode::runge_kutta4<PState>>);

  // plant_step : PState → PState
  const auto plant_step = do_step(plant, stepper, dt);

  // ---

  auto [sPlantState_unlisten, plantRecord] = util::make_listener(sPlantState);

  {
    PState x = {0., 0., 0.};
    sPlantState.send({now, x});
    for (int k = 1; k < simTime / dt; ++k) {
      x[2] = cControl.sample().ctrlVal;
      stepper.do_step(plant, x, 0, dt);
      sPlantState.send({now + k * dts, x});
    }
  }

  sPlantState_unlisten();

  {
    constexpr double margin = 0.1;
    constexpr auto sys = [](double t) {
      return -0.86502 * std::exp(-3.9537 * t) * std::sin(4.2215 * t)
             - 0.83773 * std::exp(-3.9537 * t) * std::cos(4.2215 * t)
             - 0.16226 * std::exp(-2.0924 * t) + 0.99999;
    };

    auto frpPositions =
        util::fmap([](auto x) { return x.value[0]; }, *plantRecord);
    auto realPositions = util::fmap(
        [&sys](auto x) { return sys(util::unchrono_sec(x.time - now)); },
        *plantRecord);

#ifdef PLOT
    Gnuplot gp;
    gp << "set title 'Test C: pid(" << boost::format("%.2f") % Kp << ", "
       << boost::format("%.2f") % Ki << ", " << boost::format("%.2f") % Kd
       << ")\n";
    gp << "plot '-' u 1:2:3 title 'Acceptable margin: analytical ±"
       << boost::format("%.3f") % margin
       << "' w filledcu fs solid fc rgb '#6699FF55', '-' u 1:2 title 'Test "
          "result' w l\n";
    auto range = util::fmap(
        [&sys](auto x) {
          const double t = util::unchrono_sec(x.time - now);
          return std::make_tuple(t, sys(t) + margin, sys(t) - margin);
        },
        *plantRecord);

    gp.send1d(range);
    gp.send1d(util::fmap(
        [](const auto& x) {
          return std::make_pair(util::unchrono_sec(x.time - now), x.value[0]);
        },
        *plantRecord));
    gp << "pause mouse key\nwhile (MOUSE_CHAR ne 'q') { pause mouse "
          "key; }\n";
#endif  // PLOT

    REQUIRE(util::compareVectors(frpPositions, realPositions, margin));
  }
}

TEST_CASE(
    "Test D (Proportional-Integral-Derivative Control)—Reproduct known "
    "(analytical) result, borrowed from\n"
    "http://ctms.engin.umich.edu/CTMS/"
    "index.php?example=Introduction&section=ControlPID\n"
    "A damped-driven harmonic oscillator under the influence of a "
    "PID-controller, both with parameters defined in the test, should produce "
    "a "
    "step-response within a margin of the analytical solution.",
    "[Test D], [P-controller]") {
  constexpr Real_t Kp = 350.;
  constexpr Real_t Ki = 300.;
  constexpr Real_t Kd = 50.;

  // Setpoint to x = 1, for step response.
  const sodium::cell<PState> setPoint({1., 0., 0.});
  const sodium::stream_sink<SignalPt<PState>> sPlantState;
  const auto sError = sPlantState.snapshot(setPoint, &to_error);

  const CState u0 = {now - dts, 0., 0., 0.};
  // const auto cControl = ctrl::control_frp(pid_algebra(Kp, Ki, Kd), sError,
  // u0); Or equivalently,
  const auto cControl = sError.accum<CState>(u0, pid_algebra(Kp, Ki, Kd));

  ode::runge_kutta4<PState> stepper;
  // do_step : Plant → Stepper → double → PState → PState
  //                             ^ dt     ^ x[k]   ^ x[k+1]
  constexpr auto do_step =
      curry<4>(sim::do_step_with<ode::runge_kutta4<PState>>);

  // plant_step : PState → PState
  const auto plant_step = do_step(plant, stepper, dt);

  // ---

  auto [sPlantState_unlisten, plantRecord] = util::make_listener(sPlantState);

  {
    PState x = {0., 0., 0.};
    sPlantState.send({now, x});
    for (int k = 1; k < simTime / dt; ++k) {
      x[2] = cControl.sample().ctrlVal;
      stepper.do_step(plant, x, 0, dt);
      sPlantState.send({now + k * dts, x});
    }
  }

  sPlantState_unlisten();

  {
    constexpr double margin = 0.2;
    constexpr auto sys = [](double t) {
      return -0.88562 * std::exp(-51.774 * t) * std::sin(54.918 * t)
             - 0.93656 * std::exp(-51.774 * t) * std::cos(54.918 * t)
             - 0.044211 * std::exp(-0.95866 * t)
             - 0.019219 * std::exp(-5.4933 * t) + 0.99999;
    };

    auto frpPositions =
        util::fmap([](auto x) { return x.value[0]; }, *plantRecord);
    auto realPositions = util::fmap(
        [&sys](auto x) { return sys(util::unchrono_sec(x.time - now)); },
        *plantRecord);

#ifdef PLOT
    Gnuplot gp;
    gp << "set title 'Test D: pid(" << boost::format("%.2f") % Kp << ", "
       << boost::format("%.2f") % Ki << ", " << boost::format("%.2f") % Kd
       << ")\n";
    gp << "plot '-' u 1:2:3 title 'Acceptable margin: analytical ±"
       << boost::format("%.3f") % margin
       << "' w filledcu fs solid fc rgb '#6699FF55', '-' u 1:2 title 'Test "
          "result' w l\n";
    auto range = util::fmap(
        [&sys](auto x) {
          const double t = util::unchrono_sec(x.time - now);
          return std::make_tuple(t, sys(t) + margin, sys(t) - margin);
        },
        *plantRecord);

    gp.send1d(range);
    gp.send1d(util::fmap(
        [](const auto& x) {
          return std::make_pair(util::unchrono_sec(x.time - now), x.value[0]);
        },
        *plantRecord));
    gp << "pause mouse key\nwhile (MOUSE_CHAR ne 'q') { pause mouse "
          "key; }\n";
#endif  // PLOT

    REQUIRE(util::compareVectors(frpPositions, realPositions, margin));
  }
}
