#pragma once

#include "../Plant.hpp"

namespace util {
  template <typename Stepper>
  auto do_step_with(sim::Plant p, Stepper stepper, double dt, sim::PState x0)
      -> sim::PState {
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
}  // namespace util
