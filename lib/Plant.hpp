#pragma once

#include <iostream>

#include <array>

//          ┌ · ┐  // Position
// PState = │ · │  // Speed
//          └ · ┘  // control variable, thrown in for Boost.odeint.
using PState = std::array<double, 3>;

namespace sim {
  struct Plant {
    const double staticForce;
    const double damping;
    const double spring;

    Plant(double force, double damp, double spring)
        : staticForce(force), damping(damp), spring(spring) {}

    void operator()(const PState& x, PState& dxdt, double /*time*/) const {
      dxdt[0] = x[1];
      dxdt[1] = -spring * x[0] - damping * x[1] - x[2] + staticForce;
      dxdt[2] = 0.;  // Control variable dynamics are external to integration.
    }
  };

  template <typename Stepper>
  auto do_step_with(Plant p, Stepper stepper, double dt, PState x0) -> PState {
    stepper.do_step(p, x0, 0., dt);  // time doesn't matter.
    return x0;
  }

  template <typename Predicate, typename F, typename A>
  auto step_while(Predicate p, F f, A x0) -> decltype(f(x0)) {
    A current = x0;
    A last;
    do {
      last = current;
      current = f(last);
    } while (p(last, current));

    return current;
  }

  double lyapunov(const PState& s, const PState& setPoint = {0, 0, 0}) {
    const auto error = setPoint[0] - s[0];
    const auto lyapunovEnergy = error * error + s[1] * s[1];
    std::cout << "[ e:" << error << ", v:" << s[1] << ", V:" << lyapunovEnergy
              << "]\n";
    return lyapunovEnergy;
  }
}  // namespace sim
