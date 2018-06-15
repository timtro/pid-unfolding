#pragma once

#include <chrono>
#include <memory>
#include <optional>

#include <boost/hana.hpp>

namespace chrono = std::chrono;
using namespace std::chrono_literals;

template <typename T = double, typename Clock = chrono::steady_clock>
struct SignalPt {
  chrono::time_point<Clock> time;
  T value;
};

template <typename Real = double, typename Clock = chrono::steady_clock>
struct PIDState {
  chrono::time_point<Clock> time;
  Real errSum;
  Real error;
  Real ctrlVal;
};

template <typename Real = double, typename Clock = chrono::steady_clock>
auto pid_algebra(Real kp, Real ki, Real kd) {
  return [kp, ki, kd](PIDState<Real, Clock> prev,
                      SignalPt<Real, Clock> errSigl) -> PIDState<Real, Clock> {
    const auto deltaT =
        chrono::duration_cast<chrono::seconds>(errSigl.time - prev.time)
            .count();
    const auto errSum = prev.errSum + (errSigl.value * deltaT);
    const auto dErr = (errSigl.value - prev.error) / deltaT;
    const auto ctrlVal = kp * errSigl.value + ki * errSum + kd * dErr;

    return {errSigl.time, errSum, errSigl.value, ctrlVal};
  };
}
