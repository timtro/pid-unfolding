#include <iostream>

#include <sodium/sodium.h>
#include <boost/hana/functional/curry.hpp>
#include <boost/numeric/odeint.hpp>
#include <catch/catch.hpp>

#include "../lib/control-frp.hpp"
#include "../lib/pid.hpp"

#include "Plant.hpp"
#include "test-util.hpp"

#ifdef PLOT
#include "gnuplot-iostream.h"
#endif  // PLOT

namespace ode = boost::numeric::odeint;
using boost::hana::curry;
using Real_t = double;
using CState = PIDState<Real_t>;

template <typename Clock = chrono::steady_clock>
inline constexpr auto double_to_duration(double t) {
  return chrono::duration_cast<typename Clock::duration>(
      chrono::duration<double>(t));
}

template <typename A>
inline constexpr auto unchrono_sec(A t) {
  chrono::duration<double> tAsDouble = t;
  return tAsDouble.count();
}

const auto now = chrono::steady_clock::now();
constexpr double dt = 0.01;  // seconds.
constexpr auto dts = double_to_duration(dt);

inline SignalPt<double> to_error(const SignalPt<PState>& a, const PState& b) {
  return {a.time, a.value[0] - b[0]};
}

// do_step : Plant → Stepper → double → PState → PState
//                             ^ dt     ^ x[k]   ^ x[k+1]
constexpr auto do_step = curry<4>(sim::do_step_with<ode::runge_kutta4<PState>>);

TEST_CASE(
    "A damped-driven harmonic oscillator controlled with a P-controller, both "
    "with parameters defined below, should produce a step-response with a "
    "specific overshoot. That overshoot was computed independently with "
    "MATLAB.") {
  constexpr Real_t m = 1.;
  constexpr Real_t b = 10. / m;
  constexpr Real_t k = 20. / m;
  constexpr Real_t F = 1. / m;

  constexpr Real_t Kp = 300;
  constexpr Real_t Ki = 0.;
  constexpr Real_t Kd = 0.;

  // Setpoint to x = 1, for step response.
  const sodium::cell<PState> setPoint({1., 0., 0.});

  const sodium::stream_sink<SignalPt<PState>> sPlantState;
  const auto sError = sPlantState.snapshot(setPoint, &to_error);

  const CState u0 = {now - dts, 0., 0., 0.};
  const sim::Plant plant(F, b, k);

  ode::runge_kutta4<PState> rk4;
  ode::runge_kutta4<PState> stepper;

  // plant_step : PState → PState
  const auto plant_step = do_step(plant, rk4, dt);

  const auto cControl = ctrl::control_frp(pid_algebra(Kp, Ki, Kd), sError, u0);
  // Or equivalently,
  // const auto cControl = sError.accum<CState>(u0, pid_algebra(1., 1., 1.));

  auto plantRecord = std::make_shared<std::vector<SignalPt<PState>>>();
  const auto sPlantState_unlisten =
      sPlantState.listen([plantRecord](auto x) { plantRecord->push_back(x); });

  // ---
  {
    PState x = {0., 0., 0.};
    sPlantState.send({now, x});
    for (int k = 1; k < 2 / dt; ++k) {
      x[2] = cControl.sample().ctrlVal;
      rk4.do_step(plant, x, 0, dt);
      sPlantState.send({now + k * dts, x});
    }
  }

  sPlantState_unlisten();

  {
    const auto positions =
        util::fmap([](auto x) { return x.value[0]; }, *plantRecord);

    // Assuming that the overshoot is equivalent to the max_element.
    REQUIRE(*std::max_element(positions.cbegin(), positions.cend())
            == Approx(1.378).epsilon(0.001));
  }

#ifdef PLOT
  Gnuplot gp;
  gp << "plot '-' u 1:2 w l\n";
  gp.send1d(util::fmap(
      [](const auto& x) {
        return std::make_pair(unchrono_sec(x.time - now), x.value[0]);
      },
      *plantRecord));
  gp << "pause mouse key\nwhile (MOUSE_CHAR ne 'q') { pause mouse key; }\n";
#endif  // PLOT
}
