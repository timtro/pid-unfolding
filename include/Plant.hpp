#pragma once

#include <array>

//          ┌ · ┐  // Position
// PState = │ · │  // Speed
//          └ · ┘  // control variable, thrown in for Boost.odeint.

namespace sim {
  using PState = std::array<double, 3>;

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

  double lyapunov(const PState& s, const PState& setPoint = {0, 0, 0}) {
    const auto error = setPoint[0] - s[0];
    const auto lyapunovEnergy = error * error + s[1] * s[1];

    return lyapunovEnergy;
  }
}  // namespace sim
