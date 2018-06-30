#pragma once

#include <array>

template <typename T>
inline int sgn(T x) {
  return (T(0) < x) - (x < T(0));
}

//          ┌ · ┐  // Position
// PState = │ · │  // Speed
//          └ · ┘  // control variable, thrown in for Boost.odeint.
using PState = std::array<double, 3>;

namespace sim {
  struct Plant {
    const double staticForce;
    const double damping;

    Plant(double force, double damp) : staticForce(force), damping(damp) {}

    void operator()(const PState& x, PState& dxdt, double /*time*/) const {
      dxdt[0] = x[1];
      dxdt[1] = staticForce - sgn(x[1]) * damping * x[1] * x[1] - x[2];
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
}  // namespace sim
