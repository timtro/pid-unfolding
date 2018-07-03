#pragma once

#include <chrono>
#include <functional>
#include <memory>
#include <optional>

#include <boost/hana.hpp>

#include "SignalPt.hpp"

namespace chrono = std::chrono;
using namespace std::chrono_literals;

template <typename Real = double, typename Clock = chrono::steady_clock>
struct PIDState {
  chrono::time_point<Clock> time;
  Real errSum;
  Real error;
  Real ctrlVal;
};

template <typename Real = double, typename Clock = chrono::steady_clock>
auto pid_algebra(Real kp, Real ki, Real kd) {
  return [kp, ki, kd](SignalPt<Real, Clock> errSigl,
                      PIDState<Real, Clock> prev) -> PIDState<Real, Clock> {
    const auto deltaT =
        chrono::duration_cast<chrono::nanoseconds>(errSigl.time - prev.time)
            .count();
    constexpr int nanosecondsPerSecond = 1E9;
    const auto errSum =
        prev.errSum + (errSigl.value * deltaT / nanosecondsPerSecond);
    const auto dErr =
        (errSigl.value - prev.error) * nanosecondsPerSecond / deltaT;
    const auto ctrlVal = kp * errSigl.value + ki * errSum + kd * dErr;

    return {errSigl.time, errSum, errSigl.value, ctrlVal};
  };
}
