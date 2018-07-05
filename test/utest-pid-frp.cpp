#include <cmath>
#include <iostream>

#include <sodium/sodium.h>
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
    "A damped-driven harmonic oscillator controlled with a P-controller, both "
    "with parameters defined below, should produce a step-response with a "
    "specific overshoot. That overshoot was computed independently with "
    "MATLAB.") {
  constexpr Real_t Kp = 300.;
  constexpr Real_t Ki = 0.;
  constexpr Real_t Kd = 0.;

  // Setpoint to x = 1, for step response.
  const sodium::cell<PState> setPoint({1., 0., 0.});
  const sodium::stream_sink<SignalPt<PState>> sPlantState;
  const auto sError = sPlantState.snapshot(setPoint, &to_error);

  const CState u0 = {now - dts, 0., 0., 0.};
  const auto cControl = ctrl::control_frp(pid_algebra(Kp, Ki, Kd), sError, u0);
  // Or equivalently,
  // const auto cControl = sError.accum<CState>(u0, pid_algebra(1., 1., 1.));

  ode::runge_kutta4<PState> stepper;
  // do_step : Plant → Stepper → double → PState → PState
  //                             ^ dt     ^ x[k]   ^ x[k+1]
  constexpr auto do_step =
      curry<4>(sim::do_step_with<ode::runge_kutta4<PState>>);

  // plant_step : PState → PState
  const auto plant_step = do_step(plant, stepper, dt);

  // ---

  auto plantRecord = std::make_shared<std::vector<SignalPt<PState>>>();
  const auto sPlantState_unlisten =
      sPlantState.listen([plantRecord](auto x) { plantRecord->push_back(x); });
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

    REQUIRE(util::compareVectors(frpPositions, realPositions));

#ifdef PLOT
    Gnuplot gp;
    gp << "plot '-' u 1:2:3 w filledcu fs solid fc rgb '#66bbbbbb', '-' u 1:2 "
          "w l\n";
    auto range = util::fmap(
        [&sys](auto x) {
          const double t = util::unchrono_sec(x.time - now);
          return std::make_tuple(t, sys(t) + 0.1, sys(t) - 0.1);
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
  }
}
