#include <boost/numeric/odeint.hpp>
#include <catch/catch.hpp>
#include <iostream>
#include "../lib/pid.hpp"
#include <array>

namespace odeint = boost::numeric::odeint;

using Real_t           = double;
using RealSignalVector = std::vector<SignalPt<Real_t>>;
using CState           = PIDState<Real_t>;
using PState           = std::array<int, 2>;

const auto now    = chrono::steady_clock::now();
const auto oneSec = 1s;

template <typename T, typename A, typename F>
auto get_each(std::vector<T, A> v, F accessor) {
  using AccessedType = std::invoke_result_t<F, decltype(v.at(0))>;
  std::vector<AccessedType> result;
  result.reserve(v.size());
  for (const auto &each : v) {
    result.push_back(accessor(each));
    std::vector<T, A> result(v.size());
  }
  return result;
}

constexpr auto list_fold =
    [](const auto f, const auto x0,
       const RealSignalVector &data) -> std::vector<CState> {
  auto accumulator = x0;
  std::vector<decltype(accumulator)> xs;
  xs.reserve(data.size());

  for (const auto &datum : data) {
    accumulator = f(accumulator, datum);
    xs.push_back(accumulator);
  };
  return xs;
};

struct Plant {
  const double static_force;
  const double damping;

  Plant(double force, double damp) : static_force(force), damping(damp) {}

  void operator()(const PState &x, PState &dxdt, double /*time*/) const {
    dxdt[0] = x[1];
    dxdt[1] = -static_force - damping * x[1] * x[1];
  }
};

constexpr double dt = 0.01;
odeint::runge_kutta4<PState> rk4;
