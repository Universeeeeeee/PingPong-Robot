// SPDX-License-Identifier: Apache-2.0
// Modified for standalone release: neutral API names and scoped declarations.
// Source provenance and exact changes: SOURCE_MANIFEST.json at package root.

#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <cstddef>
#include <limits>
#include <string>

namespace pingpong_planner {

using Vec3 = Eigen::Vector3d;

struct BallPhysics {
  double drag_k = 0.1261;
  double magnus_k = 0.00444;
  double restitution_h = 0.64;
  double restitution_v = 0.9215;
  // Nakashima/Ace Coulomb prior (not venue-fit).
  double coulomb_friction_mu = 0.25;
  // Retained table-contact grip law from ball_physics_venue.yaml.  This is a
  // separate candidate from the Nakashima rolling-cap law above.
  double table_tangential_gain = 0.369;
  double table_friction_cap_mu = 2.0;
  Vec3 gravity{0.0, 0.0, -9.81};
  double radius = 0.02;
};

struct TableParams {
  double length = 2.74;
  double width = 1.525;
  double y_max = 0.0;
  double net_x = 1.37;
  double net_height = 0.1525;
  double net_overhang = 0.15;
};

struct PlannerConfig {
  double x_hit = 0.15;
  Vec3 target_land{2.055, -0.7625, 0.0};
  double delta_t_flight = 0.50;
  double restitution_racket = 0.654;
  double restitution_exp_g1 = 0.759;
  double restitution_exp_g2 = -0.0441;
  double integrate_dt_s = 0.001;
  double max_predict_time_s = 2.0;
  bool adaptive_predict_horizon = false;
  double max_predict_time_cap_s = 3.0;
};

enum class SpinPhysicsMode {
  kLegacyNoSpin,
  kCoulombBounce,
  kCoulombBounceAndMagnus,
  kVenueGripBounce,
  kVenueGripBounceAndMagnus,
};

struct BallState {
  Vec3 position = Vec3::Zero();
  Vec3 velocity = Vec3::Zero();
  double source_time_s = 0.0;
  std::size_t sample_count = 0;
  double sample_span_s = 0.0;
  double residual_rms_m = std::numeric_limits<double>::quiet_NaN();
  double residual_max_m = std::numeric_limits<double>::quiet_NaN();
  // Bounce-transition fields are estimator audit only. They are not packed
  // into schema-2 and never participate in command admission.
  bool bounce_transition_used = false;
  bool bounce_epoch_active = false;
  double bounce_source_time_s = std::numeric_limits<double>::quiet_NaN();
  std::size_t pre_bounce_samples = 0;
  std::size_t post_bounce_samples = 0;
  bool valid = false;
  std::string reason = "estimator_not_ready";
};

struct StrikeTarget {
  Vec3 ball_position = Vec3::Zero();
  Vec3 ball_velocity = Vec3::Zero();
  double strike_source_time_s = 0.0;
  int predicted_bounces = 0;
  bool valid = false;
  std::string reason = "not_run";
};

struct RacketCommand {
  Vec3 position = Vec3::Zero();
  Vec3 velocity = Vec3::Zero();
  Vec3 normal{1.0, 0.0, 0.0};
  Vec3 outgoing_ball_velocity = Vec3::Zero();
  Vec3 target_land = Vec3::Zero();
  double strike_source_time_s = 0.0;
  int predicted_bounces = 0;
  bool clears_net = false;
  bool bypasses_net_posts = false;
  bool valid = false;
  std::string reason = "not_run";
};

}  // namespace pingpong_planner
