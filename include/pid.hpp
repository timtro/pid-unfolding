#pragma once

#include <boost/hana.hpp>
#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <optional>

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
    const chrono::duration<double> deltaT = errSigl.time - prev.time;
    if (deltaT <= chrono::seconds{0}) return prev;
    const auto errSum = std::fma(errSigl.value, deltaT.count(), prev.errSum);
    const auto dErr = (errSigl.value - prev.error) / deltaT.count();
    const auto ctrlVal = kp * errSigl.value + ki * errSum + kd * dErr;

    return {errSigl.time, errSum, errSigl.value, ctrlVal};
  };
}
