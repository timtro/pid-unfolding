#pragma once

#include <cmath>

namespace analyt {

  double test_A(double t) {
    return std::exp(-5. * t)
               * (-(150. * std::sin(std::sqrt(285.) * t))
                      / (31. * std::sqrt(285.))
                  - (30. * std::cos(std::sqrt(285.) * t)) / 31.)
           + 30. / 31.;
  }

  double test_B(double t) {
    return 0.10118 * std::exp(-11.145 * t) * std::sin(15.511 * t)
           - 1.0944 * std::exp(-11.145 * t) * std::cos(15.511 * t)
           + 0.15697 * std::exp(-87.708 * t) + 0.93749;
  }

  double test_C(double t) {
    return -0.86502 * std::exp(-3.9537 * t) * std::sin(4.2215 * t)
           - 0.83773 * std::exp(-3.9537 * t) * std::cos(4.2215 * t)
           - 0.16226 * std::exp(-2.0924 * t) + 0.99999;
  }

  double test_D(double t) {
    return -0.88562 * std::exp(-51.774 * t) * std::sin(54.918 * t)
           - 0.93656 * std::exp(-51.774 * t) * std::cos(54.918 * t)
           - 0.044211 * std::exp(-0.95866 * t)
           - 0.019219 * std::exp(-5.4933 * t) + 0.99999;
  }

}  // namespace analyt
