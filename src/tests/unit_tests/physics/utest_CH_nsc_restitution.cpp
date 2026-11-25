// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2026 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
//
// Test for the Newton restitution model of NSC contacts.
//
// With a coefficient of restitution e, the normal velocity right after an impact must be
// e times the impact velocity: v_rebound = e * v_impact.
//
// - drop: a sphere is dropped from rest from a height h onto a fixed box, so that
//   v_impact = sqrt(2 * g * h).
// - approach: without gravity, a sphere approaches a fixed box with speed v0 and the contact
//   closes within the first step, at a gap larger than the distance travelled at the rebound
//   speed in one step (e * v0 * dt < gap < v0 * dt). The rebound speed is e * v0.
//
// =============================================================================

#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/collision/ChCollisionModel.h"

#include "gtest/gtest.h"

using namespace chrono;

class NSCRestitutionTest : public ::testing::TestWithParam<double> {};

TEST_P(NSCRestitutionTest, drop) {
    double e = GetParam();  // coefficient of restitution
    double g = 9.81;        // gravitational acceleration
    double h = 0.5;         // drop height (gap between sphere and box at release)
    double radius = 0.1;    // sphere radius

    ChCollisionModel::SetDefaultSuggestedEnvelope(0.002);
    ChCollisionModel::SetDefaultSuggestedMargin(0.001);

    ChSystemNSC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, -g));
    sys.SetMinBounceSpeed(0.01);

    auto mat = chrono_types::make_shared<ChContactMaterialNSC>();
    mat->SetFriction(0);
    mat->SetRestitution((float)e);

    auto ground = chrono_types::make_shared<ChBodyEasyBox>(2, 2, 0.2, 1000, false, true, mat);
    ground->SetPos(ChVector3d(0, 0, -0.1));
    ground->SetFixed(true);
    sys.Add(ground);

    auto ball = chrono_types::make_shared<ChBodyEasySphere>(radius, 1000, false, true, mat);
    ball->SetPos(ChVector3d(0, 0, radius + h));
    sys.Add(ball);

    double time_step = 1e-4;
    double v_impact = 0;
    double v_rebound = 0;

    // Integrate until the sphere moves upwards after the first impact
    while (sys.GetChTime() < 1) {
        double vz_old = ball->GetPosDt().z();
        sys.DoStepDynamics(time_step);
        double vz = ball->GetPosDt().z();
        if (vz_old < 0 && vz > 0) {
            v_impact = -vz_old;
            v_rebound = vz;
            break;
        }
    }

    // Impact velocity from free fall (the contact envelope shortens the fall by a small amount)
    ASSERT_GT(v_impact, 0);
    EXPECT_NEAR(v_impact, std::sqrt(2 * g * h), 1e-2);

    // Rebound velocity from the Newton restitution model
    EXPECT_NEAR(v_rebound / v_impact, e, 1e-2);
}

TEST_P(NSCRestitutionTest, approach) {
    double e = GetParam();  // coefficient of restitution
    double v0 = 1;          // approach speed
    double radius = 0.1;    // sphere radius
    double time_step = 1e-3;
    double gap = 0.5 * (1 + e) * v0 * time_step;  // initial gap, with e * v0 * dt < gap < v0 * dt

    ChCollisionModel::SetDefaultSuggestedEnvelope(0.002);
    ChCollisionModel::SetDefaultSuggestedMargin(0.001);

    ChSystemNSC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    sys.SetGravitationalAcceleration(VNULL);
    sys.SetMinBounceSpeed(0.01);

    auto mat = chrono_types::make_shared<ChContactMaterialNSC>();
    mat->SetFriction(0);
    mat->SetRestitution((float)e);

    auto ground = chrono_types::make_shared<ChBodyEasyBox>(2, 2, 0.2, 1000, false, true, mat);
    ground->SetPos(ChVector3d(0, 0, -0.1));
    ground->SetFixed(true);
    sys.Add(ground);

    auto ball = chrono_types::make_shared<ChBodyEasySphere>(radius, 1000, false, true, mat);
    ball->SetPos(ChVector3d(0, 0, radius + gap));
    ball->SetPosDt(ChVector3d(0, 0, -v0));
    sys.Add(ball);

    sys.DoStepDynamics(time_step);

    EXPECT_NEAR(ball->GetPosDt().z(), e * v0, 1e-4);
}

INSTANTIATE_TEST_SUITE_P(Physics, NSCRestitutionTest, ::testing::Values(0.3, 0.5, 0.8));
