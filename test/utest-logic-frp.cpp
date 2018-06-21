#include <sodium/sodium.h>
#include <catch/catch.hpp>
#include <iostream>

#include "../lib/pid.hpp"
#include "Plant.hpp"

using Real_t = double;
using CState = PIDState<Real_t>;

template <typename Clock = chrono::steady_clock>
inline constexpr auto double_to_duration(double t) {
  return chrono::duration_cast<typename Clock::duration>(
      chrono::duration<double>(t));
}

constexpr double dt = 1;  // seconds.
constexpr auto dts = double_to_duration(dt);

constexpr auto add_alg = [](SignalPt<Real_t> errSigl,
                            PIDState<Real_t> prev) -> PIDState<Real_t> {
  const auto newTime = prev.time + (errSigl.time - prev.time);
  return {newTime, errSigl.value + prev.errSum + 1, prev.error + 1,
          prev.ctrlVal + 1};
};

template <typename Alg>
struct PidFrp {
  const Alg alg; // alg : (SignalPt<A>, PIDState<A>) → PIDState<A>, A = Real_t
  const sodium::stream<PState> plantOutput;
  const chrono::time_point<chrono::steady_clock> now =
      chrono::steady_clock::now();
  std::shared_ptr<std::vector<SignalPt<double>>> errorRecord{
      new std::vector<SignalPt<double>>()};
  std::shared_ptr<std::vector<CState>> ctrlScan{new std::vector<CState>()};

  PidFrp(Alg alg) : alg(alg) {}

  void operator()() {
    CState u0 = {now, 0., 0., 0.};
    sodium::cell<double> setpoint(0.);

    const auto error = plantOutput.snapshot(
        setpoint, [](const SignalPt<PState>& x, const double& setpt) {
          return SignalPt<double>{x.time, setpt - x.value[0]};
        });

    auto error_unlisten =
        error.listen([errorRecord(this->errorRecord)](auto x) {
          errorRecord->push_back(x);
        });

    const auto control = error.template accum<CState>(u0, alg);
    auto ctrl_unlisten = control.listen(
        [ctrlScan(this->ctrlScan)](auto x) { ctrlScan->push_back(x); });
  }
};

TEST_CASE("") {}
