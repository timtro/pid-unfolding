#include <sodium/sodium.h>
#include <boost/hana/functional/curry.hpp>
#include <boost/numeric/odeint.hpp>
#include <catch/catch.hpp>

#include "../include/Plant.hpp"

using boost::hana::curry;
using sim::PState;
namespace odeint = boost::numeric::odeint;

TEST_CASE(
    "Study terminal velocity without spring term. Starting from rest at the "
    "origin …") {
  // The plant models a point mass with a static force and a quadratic drag
  // term.
  // (*)  Fs - η v = m dv/dt
  //     where
  //       * Fs is the static force,
  //       * η is the drag coefficient,
  //       * m is the inertial mass, and
  //       * v is the linear speed.
  // The terminal velocity is apparent in (*), since dv/dt will be 0:
  // (†) v = (Fs / η)

  constexpr auto do_step =
      curry<4>(sim::do_step_with<odeint::runge_kutta4<PState>>);
  constexpr double dt = 0.1;
  constexpr PState x0{0., 0., 0.};

  constexpr double convergenceTolerance = 1E-12;
  constexpr auto v_diff_not_within_tolerance = [](PState last, PState current) {
    return (current[1] - last[1]) > convergenceTolerance;
  };

  SECTION("and given (Fs, η) = (4, 1), the terminal velocity should be ≈ 4.") {
    const auto plant = sim::Plant(4., 1., 0.);
    odeint::runge_kutta4<PState> stepper;

    const auto plant_step = do_step(plant, stepper, dt);
    PState atTerminalVel =
        step_while(v_diff_not_within_tolerance, plant_step, x0);
    REQUIRE(atTerminalVel[1] == Approx(4));
  }

  SECTION("and given (Fs, η) = (9, 1), the terminal velocity should be ≈ 9.") {
    const auto plant = sim::Plant(9., 1., 0.);
    odeint::runge_kutta4<PState> stepper;

    const auto plant_step = do_step(plant, stepper, dt);
    PState atTerminalVel =
        step_while(v_diff_not_within_tolerance, plant_step, x0);
    REQUIRE(atTerminalVel[1] == Approx(9));
  }
}
