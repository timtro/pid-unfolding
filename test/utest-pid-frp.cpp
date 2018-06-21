#include <sodium/sodium.h>
#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>
#include <catch/catch.hpp>
#include <iostream>
#include "../lib/pid.hpp"

#include "Plant.hpp"

namespace odeint = boost::numeric::odeint;
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

odeint::runge_kutta4<PState> rk4;

TEST_CASE("") {
  Plant p{0., 0.};
  CState u0 = {now - dts, 0., 0., 0.};

  // plantOutput are just perfect state measurements.
  sodium::stream_sink<SignalPt<PState>> plantOutput;
  sodium::cell<double> setpoint(0.);

  const auto error = plantOutput.snapshot(
      setpoint, [](const SignalPt<PState>& x, const double& setpt) {
        return SignalPt<double>{x.time, setpt - x.value[0]};
      });
  auto errorRecord = std::make_shared<std::vector<SignalPt<double>>>();
  auto error_unlisten =
      error.listen([errorRecord](auto x) { errorRecord->push_back(x); });

  const auto control = error.accum<CState>(u0, pid_algebra(1., 1., 1.));
  auto ctrlScan = std::make_shared<std::vector<CState>>();
  auto ctrl_unlisten =
      control.listen([ctrlScan](auto x) { ctrlScan->push_back(x); });

  std::vector<PState> simulation;

  {
    PState x = {0., 0., 0.};
    simulation.push_back(x);
    plantOutput.send({now, x});
    for (double t = dt; t < 1; t += dt) {
      x[2] = control.sample().ctrlVal;
      rk4.do_step(p, x, t, dt);
      simulation.push_back(x);
      plantOutput.send({now + double_to_duration(t), x});
    }
  }

  ctrl_unlisten();
  error_unlisten();

  // for (auto& each : simulation)
  //   std::cout << each[0] << ", " << each[1] << ", " << each[2] << std::endl;

  // for (auto& each : *ctrlScan)
  //   std::cout << "Time: " << each.time.time_since_epoch().count()
  //             << " Error: " << each.error
  //             << " errSum: " << each.errSum
  //             << " ctrlVal: " << each.ctrlVal << std::endl;

  // for (auto& each: *errorRecord)
  //   std::cout << "Time: " << each.time.time_since_epoch().count()
  //             << ", Error: " << each.value
  //             << std::endl;
}
