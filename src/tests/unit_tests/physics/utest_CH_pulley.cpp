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
// Test for ChLinkLockPulley.
//
// Two wheels with radii r1 and r2, on revolute joints to ground, are connected by a belt.
// Wheel 1 is driven at constant angular speed w1 and a constant torque T is applied to wheel 2.
// For an inextensible belt, the analytical solution is:
//    w2    = w1 * r1 / r2
//    Drot2 = Drot1 * r1 / r2
//    belt tension  F = T / r2
//    motor torque |M| = F * r1 = T * r1 / r2
//
// =============================================================================

#include "chrono/physics/ChLinkLockPulley.h"
#include "chrono/physics/ChLinkMotorRotationSpeed.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/solver/ChDirectSolverLS.h"

#include "gtest/gtest.h"

using namespace chrono;

TEST(ChLinkLockPulley, transmission) {
    double r1 = 1.0;  // radius of wheel 1
    double r2 = 1.5;  // radius of wheel 2
    double w1 = 6;    // angular speed of wheel 1
    double T = 10;    // torque applied to wheel 2

    ChSystemNSC sys;
    sys.SetGravitationalAcceleration(VNULL);
    sys.SetSolver(chrono_types::make_shared<ChSolverSparseQR>());

    auto ground = chrono_types::make_shared<ChBody>();
    ground->SetFixed(true);
    sys.Add(ground);

    auto wheel1 = chrono_types::make_shared<ChBody>();
    wheel1->SetMass(10);
    wheel1->SetInertiaXX(ChVector3d(1, 1, 1));
    wheel1->SetPos(ChVector3d(0, 0, 0));
    sys.Add(wheel1);

    auto wheel2 = chrono_types::make_shared<ChBody>();
    wheel2->SetMass(10);
    wheel2->SetInertiaXX(ChVector3d(1, 1, 1));
    wheel2->SetPos(ChVector3d(0, -5, 0));
    sys.Add(wheel2);
    unsigned int acc2 = wheel2->AddAccumulator();
    wheel2->AccumulateTorque(acc2, ChVector3d(0, 0, T), false);  // constant torque, never cleared

    auto motor = chrono_types::make_shared<ChLinkMotorRotationSpeed>();
    motor->Initialize(wheel1, ground, ChFramed(wheel1->GetPos()));
    motor->SetSpeedFunction(chrono_types::make_shared<ChFunctionConst>(w1));
    sys.AddLink(motor);

    auto revolute2 = chrono_types::make_shared<ChLinkLockRevolute>();
    revolute2->Initialize(wheel2, ground, ChFramed(wheel2->GetPos()));
    sys.AddLink(revolute2);

    auto pulley = chrono_types::make_shared<ChLinkLockPulley>();
    pulley->Initialize(wheel1, wheel2, ChFramed());
    pulley->SetRadius1(r1);
    pulley->SetRadius2(r2);
    sys.AddLink(pulley);

    // Place the pulley markers (set at updates of visual assets) before the first step
    sys.Setup();
    sys.Update();

    // total rotations are measured from the interaxis direction
    double rot1_0 = pulley->GetRotation1();
    double rot2_0 = pulley->GetRotation2();

    double time_step = 1e-3;
    while (sys.GetChTime() < 2) {
        sys.DoStepDynamics(time_step);

        // the belt length (imposed distance between the belt ends) must be maintained
        ASSERT_LT(pulley->GetConstraintViolation().lpNorm<Eigen::Infinity>(), 1e-6);
    }

    // speed and rotation ratio
    double w1_sim = wheel1->GetAngVelParent().z();
    double w2_sim = wheel2->GetAngVelParent().z();
    EXPECT_NEAR(w1_sim, w1, 1e-6);
    EXPECT_NEAR(w2_sim, w1 * r1 / r2, 1e-4);
    EXPECT_NEAR(pulley->GetRotation2() - rot2_0, (pulley->GetRotation1() - rot1_0) * r1 / r2, 1e-3);

    // motor torque balancing the torque on wheel 2 through the belt
    EXPECT_NEAR(std::abs(motor->GetMotorTorque()), T * r1 / r2, 1e-3);
}
