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
// Test for SMC contacts with stiff contact Jacobians (ChSystemSMC::SetContactStiff).
//
// A sphere of mass m rests on a fixed box. The contact is a linear spring-dashpot with user
// coefficients kn and gn, with a time step h such that sqrt(kn / m) * h > 2, for which the
// explicitly integrated contact force is unstable. With the contact Jacobians included in the
// implicit integration, the sphere must settle at the static penetration m * g / kn.
//
// =============================================================================

#include <cmath>

#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/solver/ChDirectSolverLS.h"

#include "gtest/gtest.h"

using namespace chrono;

TEST(SMCStiffContact, resting_sphere) {
    double g = 9.81;
    double radius = 0.1;
    double kn = 1e7;
    double h = 2e-3;

    ChSystemSMC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, -g));
    sys.UseMaterialProperties(false);
    sys.SetContactForceModel(ChSystemSMC::Hooke);
    sys.SetContactStiff(true);
    sys.SetSolver(chrono_types::make_shared<ChSolverSparseLU>());
    sys.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);

    auto ball_mat = chrono_types::make_shared<ChContactMaterialSMC>();
    auto ground_mat = chrono_types::make_shared<ChContactMaterialSMC>();

    auto ground = chrono_types::make_shared<ChBodyEasyBox>(2, 2, 0.2, 1000, false, true, ground_mat);
    ground->SetPos(ChVector3d(0, 0, -0.1));
    ground->SetFixed(true);
    sys.Add(ground);

    auto ball = chrono_types::make_shared<ChBodyEasySphere>(radius, 1000, false, true, ball_mat);
    ball->SetPos(ChVector3d(0, 0, radius + 0.01));
    sys.Add(ball);

    // Linear contact with damping ratio 0.5 (the damping coefficient is given per unit mass)
    double m = ball->GetMass();
    for (auto mat : {ball_mat, ground_mat}) {
        mat->SetFriction(0.5f);
        mat->SetKn((float)kn);
        mat->SetGn((float)(2 * 0.5 * std::sqrt(kn / m)));
        mat->SetKt((float)kn);
        mat->SetGt(0);
    }
    ASSERT_GT(std::sqrt(kn / m) * h, 2.0);

    double z_static = radius - m * g / kn;
    double max_dev = 0;
    while (sys.GetChTime() < 1.0) {
        sys.DoStepDynamics(h);
        if (sys.GetChTime() > 0.5)
            max_dev = std::max(max_dev, std::abs(ball->GetPos().z() - z_static));
    }
    std::cout << "static penetration " << m * g / kn << "  max deviation " << max_dev << std::endl;
    ASSERT_LT(max_dev, 1e-6);
}
