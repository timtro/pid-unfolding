#pragma once

#include <sodium/sodium.h>

#include "SignalPt.hpp"

namespace ctrl {

  template <typename Alg, typename X, typename U>
  [[nodiscard]] inline auto control_scan(Alg alg,
                                         const sodium::stream<X>& cError,
                                         const U u0) {
    return cError.template accum<U>(u0, alg);
  }

  template <typename X, typename U, typename Alg, typename Fdiff>
  [[nodiscard]] inline auto control_frp(
      Alg alg, Fdiff differencer, const sodium::cell<X>& setPoint,
      const sodium::stream<SignalPt<X>>& sPlantState, U u0) {
    const auto sError = sPlantState.snapshot(setPoint, differencer);
    const auto cControl = sError.template accum<U>(u0, alg);

    return cControl;
  }

}  // namespace ctrl
