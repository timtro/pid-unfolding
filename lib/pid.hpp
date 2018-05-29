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
  Real setPt;
  Real errSum;
  Real error;
  Real ctrlVal;
};

template <typename Real = double, typename Clock = chrono::steady_clock>
auto pid_algebra(Real kp, Real ki, Real kd) {
  return [kp, ki, kd](PIDState<Real, Clock> s, SignalPt<Real, Clock> input)
    -> PIDState<Real, Clock> {
    const auto delta = s.time - input.time;
    const auto err = s.setPt - input.value;
    const auto errSum =
        s.errSum
        + (err * chrono::duration_cast<chrono::seconds>(delta).count());
    const auto dErr =
        (err - s.error) / chrono::duration_cast<chrono::seconds>(delta).count();
    const auto ctrlVal = kp * err + ki * errSum + kd * dErr;

    return {input.time, s.setPt, errSum, err, ctrlVal};
  };
}
