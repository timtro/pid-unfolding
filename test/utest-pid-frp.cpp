#include <iostream>

#include <sodium/sodium.h>
#include <boost/hana/functional/curry.hpp>
#include <boost/numeric/odeint.hpp>
#include <catch/catch.hpp>

#include "../lib/control-frp.hpp"
#include "../lib/pid.hpp"

#include "Plant.hpp"
#include "test-util.hpp"

#include "gnuplot-iostream.h"

namespace ode = boost::numeric::odeint;
using boost::hana::curry;
using Real_t = double;
using CState = PIDState<Real_t>;

template <typename Clock = chrono::steady_clock>
inline constexpr auto double_to_duration(double t) {
  return chrono::duration_cast<typename Clock::duration>(
      chrono::duration<double>(t));
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

TEST_CASE("") {
  const sodium::cell<PState> setPoint({1., 0., 0.});
  const sodium::stream_sink<SignalPt<PState>> sPlantState;
  const auto sError = sPlantState.snapshot(setPoint, &to_error);

  ode::runge_kutta4<PState> rk4;
  ode::runge_kutta4<PState> stepper;

  const CState u0 = {now - dts, 0., 0., 0.};
  const sim::Plant plant{-1., 20.};

  // plant_step : PState → PState
  const auto plant_step = do_step(plant, rk4, dt);

  const auto cControl = ctrl::control_frp(pid_algebra(1., 1., 1.), sError, u0);
  // Or equivalently,
  // const auto cControl = sError.accum<CState>(u0, pid_algebra(1., 1., 1.));

  // Listeners, to monitor and test output.
  auto errorRecord = std::make_shared<std::vector<SignalPt<Real_t>>>();
  const auto sError_unlisten = sError.listen([errorRecord](auto x) {
    errorRecord->push_back(x);
    // std::cout << "     sError(" << errorRecord->size() << "): "
    //          << " [ " << (x.time - now).count() << " ] {" << x.value << "\n";
  });

  auto ctrlRecord = std::make_shared<std::vector<CState>>();
  const auto cControl_unlisten = cControl.listen([ctrlRecord](auto x) {
    ctrlRecord->push_back(x);
    // std::cout << "   cControl(" << ctrlRecord->size() << "): "
    //           << " [ " << (x.time - now).count() << " ] "
    //           << "{errSum:" << x.errSum << ", error:" << x.error
    //           << ", ctrlVal:" << x.ctrlVal << "}\n";
  });

  auto plantRecord = std::make_shared<std::vector<SignalPt<PState>>>();
  const auto sPlantState_unlisten = sPlantState.listen([plantRecord](auto x) {
    plantRecord->push_back(x);
    // std::cout << "sPlantState(" << plantRecord->size() << "): "
    //           << " [ " << (x.time - now).count() << " ] "
    //           << "{ x:" << x.value[0] << ", v:" << x.value[1]
    //           << ", u:" << x.value[2] << "}\n";
  });

  {
    PState x = {0., 0., 0.};
    sPlantState.send({now, x});
    for (int k = 1; k < 10001; ++k) {
      x[2] = cControl.sample().ctrlVal;
      rk4.do_step(plant, x, 0, dt);
      sPlantState.send({now + k * dts, x});
    }
  }

  sError_unlisten();
  cControl_unlisten();
  sPlantState_unlisten();

  Gnuplot gp;
  gp << "plot '-' u 1:2 w l\n";
  gp.send1d(util::fmap(
      [](auto x) { return std::make_pair((x.time-now).count(), x.value[0]); },
      *plantRecord));
}
