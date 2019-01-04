#pragma once

#include <iostream>

#include <chrono>
#include <functional>
#include <memory>
#include <optional>

#include <boost/hana.hpp>

#include "SignalPt.hpp"

namespace chrono = std::chrono;
using namespace std::chrono_literals;

template <typename Clock = chrono::steady_clock>
struct PIDState {
  chrono::time_point<Clock> time;
  double errSum;
  double error;
  double ctrlVal;
};

template <typename Clock = chrono::steady_clock>
auto pid_algebra(double kp, double ki, double kd) {
  return [kp, ki, kd](PIDState<Clock> prev,
                      SignalPt<double, Clock> errSigl) -> PIDState<Clock> {
    const auto deltaT =
        chrono::duration_cast<chrono::nanoseconds>(errSigl.time - prev.time)
            .count();
    if (deltaT <= 0)
      return prev;
    constexpr int nanosecondsPerSecond = 1E9;
    const auto errSum =
        prev.errSum + (errSigl.value * deltaT / nanosecondsPerSecond);
    const auto dErr =
        (errSigl.value - prev.error) * nanosecondsPerSecond / deltaT;
    const auto ctrlVal = kp * errSigl.value + ki * errSum + kd * dErr;

    return {errSigl.time, errSum, errSigl.value, ctrlVal};
  };
}
