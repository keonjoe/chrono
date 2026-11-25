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
// Test for the APGD solver on a hanging chain.
//
// A chain of N bodies of mass m, connected by spherical joints and suspended from ground,
// hangs at rest under gravity. Joint i (counted from the top, i = 0..N-1) carries the
// weight of the N-i bodies below it:
//    |F_i| = (N - i) * m * g
// and the chain must remain at rest.
//
// =============================================================================

#include "chrono/physics/ChLinkLock.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/solver/ChSolverAPGD.h"

#include "gtest/gtest.h"

using namespace chrono;

TEST(ChSolverAPGD, hanging_chain) {
    int N = 10;       // number of bodies
    double m = 1;     // body mass
    double g = 9.81;  // gravitational acceleration
    double L = 1;     // body length

    ChSystemNSC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, -g, 0));

    auto solver = chrono_types::make_shared<ChSolverAPGD>();
    solver->SetMaxIterations(500);
    solver->SetTolerance(1e-10);
    solver->EnableWarmStart(false);
    sys.SetSolver(solver);

    auto ground = chrono_types::make_shared<ChBody>();
    ground->SetFixed(true);
    sys.Add(ground);

    std::vector<std::shared_ptr<ChLinkLockSpherical>> joints;
    std::vector<std::shared_ptr<ChBody>> bodies;
    auto prev = ground;
    for (int i = 0; i < N; i++) {
        auto body = chrono_types::make_shared<ChBody>();
        body->SetMass(m);
        body->SetInertiaXX(ChVector3d(0.1, 0.1, 0.1));
        body->SetPos(ChVector3d(0, -(i + 0.5) * L, 0));
        sys.Add(body);

        auto joint = chrono_types::make_shared<ChLinkLockSpherical>();
        joint->Initialize(body, prev, ChFramed(ChVector3d(0, -i * L, 0)));
        sys.AddLink(joint);

        bodies.push_back(body);
        joints.push_back(joint);
        prev = body;
    }

    double time_step = 1e-3;
    for (int k = 0; k < 10; k++)
        sys.DoStepDynamics(time_step);

    for (int i = 0; i < N; i++) {
        double F = joints[i]->GetReaction2().force.Length();
        EXPECT_NEAR(F, (N - i) * m * g, 1e-3 * N * m * g) << "joint " << i;
        EXPECT_NEAR(bodies[i]->GetPosDt().Length(), 0, 1e-6) << "body " << i;
    }
}
