#pragma once

#include <sodium/sodium.h>

namespace ctrl {

  template <typename Alg, typename X, typename U>
  [[nodiscard]] inline auto control_frp(Alg alg,
                                        const sodium::stream<X>& cError,
                                        const U u0) {
    return cError.template accum<U>(u0, alg);
  }

}  // namespace ctrl
