#pragma once

#include <cmath>

using std::sqrt;
using std::sin;
using std::cos;
using std::exp;

namespace analyt {

  double test_A(double t) {
    return -0.27291*exp(-5*t)*sin(17.175*t)-0.9375*exp(-5*t)*cos(17.175*t)+0.9375;
  }

  double test_B(double t) {
    return 0.042137*exp(-10*t)*sin(14.832*t)-0.9375*exp(-10*t)*cos(14.832*t)+0.9375;
  }

  double test_C(double t) {
    return -0.86502*exp(-3.9537*t)*sin(4.2215*t)-0.83773*exp(-3.9537*t)*cos(4.2215*t)-0.16226*exp(-2.0924*t)+0.99999;
  }

  double test_D(double t) {
    return -0.043992*exp(-0.95693*t)-0.017952*exp(-5.899*t)-0.93805*exp(-53.144*t)+1.0;
  }

}  // namespace analyt
