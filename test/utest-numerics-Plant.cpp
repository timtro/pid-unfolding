#include <sodium/sodium.h>
#include <boost/hana/functional/curry.hpp>
#include <boost/numeric/odeint.hpp>
#include <catch/catch.hpp>

#include "Plant.hpp"

using boost::hana::curry;
namespace odeint = boost::numeric::odeint;

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

TEST_CASE("Study terminal velocity. Starting from rest at the origin …") {
  // The plant models a point mass with a static force and a quadratic drag
  // term.
  // (*)  Fs - η v^2 = m dv/dt
  //     where
  //       * Fs is the static force,
  //       * η is the drag coefficient,
  //       * m is the inertial mass, and
  //       * v is the linear speed.
  // The terminal velocity is apparent in (*), since dv/dt will be 0:
  // (†) v = (Fs / η)^(1/2)

  constexpr auto do_step = curry<4>(do_step_with<odeint::runge_kutta4<PState>>);
  constexpr double dt = 0.1;
  constexpr PState x0{0., 0., 0.};

  constexpr double convergenceTolerance = 1E-12;
  constexpr auto v_diff_not_within_tolerance = [](PState last, PState current) {
    return (current[1] - last[1]) > convergenceTolerance;
  };

  SECTION(
      "and given (Fs, η) = (4, 1), the terminal velocity should be ≈ 2.") {
    const auto plant = Plant(4., 1.);
    odeint::runge_kutta4<PState> stepper;

    const auto step = do_step(plant, stepper, dt);
    PState atTerminalVel = step_while(v_diff_not_within_tolerance, step, x0);
    REQUIRE(atTerminalVel[1] == Approx(2));
  }

  SECTION(
      "and given (Fs, η) = (9, 1), the terminal velocity should be ≈ 3.") {
    const auto plant = Plant(9., 1.);
    odeint::runge_kutta4<PState> stepper;

    const auto step = do_step(plant, stepper, dt);
    PState atTerminalVel = step_while(v_diff_not_within_tolerance, step, x0);
    REQUIRE(atTerminalVel[1] == Approx(3));
  }
}
