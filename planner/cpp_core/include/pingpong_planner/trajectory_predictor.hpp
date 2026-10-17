// SPDX-License-Identifier: Apache-2.0
// Modified for standalone release: neutral API names and scoped declarations.
// Source provenance and exact changes: SOURCE_MANIFEST.json at package root.

#pragma once

#include "pingpong_planner/types.hpp"

namespace pingpong_planner {

class TrajectoryPredictor {
 public:
  TrajectoryPredictor(BallPhysics physics, PlannerConfig config, TableParams table);

  StrikeTarget predict(const BallState& state, double x_hit) const noexcept;

  // Spin-aware research path for offline prediction.
  // kLegacyNoSpin dispatches to predict() exactly.  Callers decide whether a
  // prediction is usable; this method does not perform command admission.
  StrikeTarget predict_with_spin(
      const BallState& state,
      double x_hit,
      const Vec3& omega_rad_s,
      SpinPhysicsMode mode) const noexcept;

 private:
  double prediction_horizon_s(const BallState& state, double x_hit) const noexcept;
  Vec3 flight_acceleration(
      const Vec3& velocity,
      const Vec3& omega_rad_s,
      bool use_magnus) const noexcept;
  void apply_coulomb_bounce(
      const Vec3& velocity_in,
      const Vec3& omega_in,
      Vec3* velocity_out,
      Vec3* omega_out) const noexcept;
  void apply_venue_grip_bounce(
      const Vec3& velocity_in,
      const Vec3& omega_in,
      Vec3* velocity_out,
      Vec3* omega_out) const noexcept;

  BallPhysics physics_;
  PlannerConfig config_;
  TableParams table_;
};

}  // namespace pingpong_planner
