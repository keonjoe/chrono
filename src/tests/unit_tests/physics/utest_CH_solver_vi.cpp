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
// Tests for the iterative VI solvers on a hanging chain.
//
// A chain of N bodies of masses m_j, connected by spherical joints and suspended
// from ground, hangs at rest under gravity. Joint i (counted from the top) carries
// the weight of the bodies below it:
//    |F_i| = g * sum_{j >= i} m_j
//
// - BB: a heavy chain with a mass ratio of 1000 between neighbors. The reactions
//   must be accurate to the solver tolerance (the projected gradient used for the
//   stopping test must not vanish by cancellation while the multipliers are large).
// - APGD: a uniform chain must reach the tolerance well before the iteration limit
//   (the line search must not stall the step size), with and without the diagonal
//   preconditioner.
// - APGD: without constraints, the reported error is zero.
// - PJacobi: with warm start, the reaction of a single pendulum at rest stays m*g.
// - BB, APGD, PSOR, PJacobi: joints between two fixed bodies and between a fixed body and ground give constraint rows
//   that act on no active variable (zero Schur complement). They must not affect the solution: a box resting on the
//   fixed floor carries m*g, a ball spun about the vertical by a rotational speed motor and standing on the floor
//   (contact on the spin axis, no slip) carries m*g and turns at the imposed speed, and the joint reactions are zero.
// - PSOR: a box resting on a floor stays at rest when the solver stops at its tolerance (the stopping test must see
//   a contact that keeps a positive reaction while separating).
//
// =============================================================================

#include "chrono/functions/ChFunctionConst.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChLinkLock.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChLinkMotorRotationSpeed.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/solver/ChSolverAPGD.h"
#include "chrono/solver/ChSolverBB.h"
#include "chrono/solver/ChSolverPJacobi.h"
#include "chrono/solver/ChSolverPSOR.h"

#include "gtest/gtest.h"

using namespace chrono;

static const double g = 9.81;

// Build a hanging chain of N bodies (masses m0, m0*ratio, m0, ...) in the given system.
static void BuildChain(ChSystemNSC& sys, int N, double m0, double ratio, std::vector<std::shared_ptr<ChLinkLockSpherical>>& joints, std::vector<std::shared_ptr<ChBody>>& bodies) {
    sys.SetGravitationalAcceleration(ChVector3d(0, -g, 0));

    auto ground = chrono_types::make_shared<ChBody>();
    ground->SetFixed(true);
    sys.Add(ground);

    auto prev = ground;
    for (int i = 0; i < N; i++) {
        double m = (i % 2) ? m0 * ratio : m0;
        auto body = chrono_types::make_shared<ChBody>();
        body->SetMass(m);
        body->SetInertiaXX(ChVector3d(0.1, 0.1, 0.1) * m);
        body->SetPos(ChVector3d(0, -(i + 0.5), 0));
        sys.Add(body);

        auto joint = chrono_types::make_shared<ChLinkLockSpherical>();
        joint->Initialize(body, prev, ChFramed(ChVector3d(0, -i, 0)));
        sys.AddLink(joint);

        bodies.push_back(body);
        joints.push_back(joint);
        prev = body;
    }
}

// Largest relative error of the joint reactions with respect to the carried weight.
static double ReactionError(const std::vector<std::shared_ptr<ChLinkLockSpherical>>& joints, const std::vector<std::shared_ptr<ChBody>>& bodies) {
    double err = 0;
    for (size_t i = 0; i < joints.size(); i++) {
        double W = 0;
        for (size_t j = i; j < bodies.size(); j++)
            W += bodies[j]->GetMass() * g;
        err = std::max(err, std::abs(joints[i]->GetReaction2().force.Length() - W) / W);
    }
    return err;
}

TEST(ChSolverBB, heavy_chain) {
    ChSystemNSC sys;
    std::vector<std::shared_ptr<ChLinkLockSpherical>> joints;
    std::vector<std::shared_ptr<ChBody>> bodies;
    BuildChain(sys, 10, 1, 1000, joints, bodies);

    auto solver = chrono_types::make_shared<ChSolverBB>();
    solver->SetMaxIterations(500);
    solver->SetTolerance(1e-10);
    solver->EnableWarmStart(false);
    sys.SetSolver(solver);

    // no solve yet
    EXPECT_EQ(solver->GetError(), 0);

    for (int k = 0; k < 10; k++)
        sys.DoStepDynamics(1e-3);

    EXPECT_LT(ReactionError(joints, bodies), 1e-8);
    for (auto& body : bodies)
        EXPECT_NEAR(body->GetPosDt().Length(), 0, 1e-9);
}

TEST(ChSolverAPGD, chain_convergence) {
    for (bool precond : {true, false}) {
        ChSystemNSC sys;
        std::vector<std::shared_ptr<ChLinkLockSpherical>> joints;
        std::vector<std::shared_ptr<ChBody>> bodies;
        BuildChain(sys, 10, 1, 1, joints, bodies);

        int max_iterations = 300;
        auto solver = chrono_types::make_shared<ChSolverAPGD>();
        solver->SetMaxIterations(max_iterations);
        solver->SetTolerance(1e-10);
        solver->EnableWarmStart(false);
        solver->EnableDiagonalPreconditioner(precond);
        sys.SetSolver(solver);

        for (int k = 0; k < 10; k++) {
            sys.DoStepDynamics(1e-3);
            EXPECT_LT(solver->GetIterations(), max_iterations) << "precond " << precond << ", step " << k;
            EXPECT_LT(solver->GetError(), 1e-10) << "precond " << precond << ", step " << k;
        }

        EXPECT_LT(ReactionError(joints, bodies), 1e-12) << "precond " << precond;
        for (auto& body : bodies)
            EXPECT_NEAR(body->GetPosDt().Length(), 0, 1e-10) << "precond " << precond;
    }
}

TEST(ChSolverAPGD, no_constraints) {
    ChSystemNSC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, -g, 0));
    auto body = chrono_types::make_shared<ChBody>();
    sys.Add(body);

    auto solver = chrono_types::make_shared<ChSolverAPGD>();
    sys.SetSolver(solver);

    double h = 1e-3;
    sys.DoStepDynamics(h);
    EXPECT_EQ(solver->GetError(), 0);
    EXPECT_NEAR(body->GetPosDt().y(), -g * h, 1e-12);
}

TEST(ChSolverPJacobi, warm_start) {
    ChSystemNSC sys;
    std::vector<std::shared_ptr<ChLinkLockSpherical>> joints;
    std::vector<std::shared_ptr<ChBody>> bodies;
    BuildChain(sys, 1, 1, 1, joints, bodies);

    auto solver = chrono_types::make_shared<ChSolverPJacobi>();
    solver->SetMaxIterations(200);
    solver->SetTolerance(1e-10);
    solver->EnableWarmStart(true);
    sys.SetSolver(solver);

    for (int k = 0; k < 10; k++) {
        sys.DoStepDynamics(1e-3);
        EXPECT_LT(ReactionError(joints, bodies), 1e-6) << "step " << k;
    }
    EXPECT_NEAR(bodies[0]->GetPosDt().Length(), 0, 1e-9);
}

TEST(ChSolverVI, links_between_fixed_bodies) {
    std::vector<std::shared_ptr<ChIterativeSolverVI>> solvers = {chrono_types::make_shared<ChSolverBB>(), chrono_types::make_shared<ChSolverAPGD>(),
                                                                 chrono_types::make_shared<ChSolverPSOR>(), chrono_types::make_shared<ChSolverPJacobi>()};
    for (auto& solver : solvers) {
        std::string name = solver->GetTypeAsString();

        ChSystemNSC sys;
        sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
        sys.SetGravitationalAcceleration(ChVector3d(0, 0, -g));
        solver->SetMaxIterations(600);
        solver->SetTolerance(1e-10);
        sys.SetSolver(solver);

        auto mat = chrono_types::make_shared<ChContactMaterialNSC>();
        mat->SetFriction(0.4f);

        auto ground = chrono_types::make_shared<ChBody>();
        ground->SetFixed(true);
        sys.Add(ground);

        auto floor = chrono_types::make_shared<ChBodyEasyBox>(4, 4, 0.2, 1000, false, true, mat);
        floor->SetPos(ChVector3d(0, 0, -0.1));
        floor->SetFixed(true);
        sys.Add(floor);

        auto wall = chrono_types::make_shared<ChBodyEasyBox>(0.2, 4, 1, 1000, false, false, mat);
        wall->SetPos(ChVector3d(2.1, 0, 0.5));
        wall->SetFixed(true);
        sys.Add(wall);

        // joints between a fixed body and ground and between two fixed bodies
        auto fix1 = chrono_types::make_shared<ChLinkMateFix>();
        fix1->Initialize(floor, ground, ChFramed(ChVector3d(0, 0, -0.1)));
        sys.Add(fix1);
        auto fix2 = chrono_types::make_shared<ChLinkMateFix>();
        fix2->Initialize(wall, floor, ChFramed(ChVector3d(2.1, 0, 0.5)));
        sys.Add(fix2);

        // box resting on the floor
        auto box = chrono_types::make_shared<ChBodyEasyBox>(0.4, 0.4, 0.4, 1000, false, true, mat);
        box->SetPos(ChVector3d(-1, 0, 0.2));
        sys.Add(box);

        // ball standing on the floor, turned by a speed motor about the vertical axis (free to move vertically)
        double w = CH_PI / 2;
        auto ball = chrono_types::make_shared<ChBodyEasySphere>(0.15, 4000, false, true, mat);
        ball->SetPos(ChVector3d(0.8, 0, 0.15));
        sys.Add(ball);
        auto motor = chrono_types::make_shared<ChLinkMotorRotationSpeed>();
        motor->Initialize(ball, ground, ChFramed(ChVector3d(0.8, 0, 0.15)));
        motor->SetSpeedFunction(chrono_types::make_shared<ChFunctionConst>(w));
        motor->SetSpindleConstraint(ChLinkMotorRotation::SpindleConstraint::CYLINDRICAL);
        sys.Add(motor);

        for (int k = 0; k < 200; k++)
            sys.DoStepDynamics(5e-3);

        double Wbox = box->GetMass() * g;
        double Wball = ball->GetMass() * g;
        EXPECT_NEAR(box->GetContactForce().z(), Wbox, 1e-3 * Wbox) << name;
        EXPECT_NEAR(ball->GetContactForce().z(), Wball, 1e-3 * Wball) << name;
        EXPECT_NEAR(box->GetPos().z(), 0.2, 1e-3) << name;
        EXPECT_NEAR(ball->GetAngVelParent().z(), w, 1e-6) << name;
        EXPECT_EQ(fix1->GetReaction2().force.Length(), 0) << name;
        EXPECT_EQ(fix2->GetReaction2().force.Length(), 0) << name;
        EXPECT_LT(solver->GetError(), 1e-3) << name;
    }
}

TEST(ChSolverPSOR, box_at_rest) {
    ChSystemNSC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, -g));

    auto solver = chrono_types::make_shared<ChSolverPSOR>();
    solver->SetMaxIterations(600);
    solver->SetTolerance(1e-10);
    sys.SetSolver(solver);

    auto mat = chrono_types::make_shared<ChContactMaterialNSC>();
    mat->SetFriction(0.4f);

    auto floor = chrono_types::make_shared<ChBodyEasyBox>(4, 4, 0.2, 1000, false, true, mat);
    floor->SetPos(ChVector3d(0, 0, -0.1));
    floor->SetFixed(true);
    sys.Add(floor);

    auto box = chrono_types::make_shared<ChBodyEasyBox>(0.4, 0.4, 0.4, 1000, false, true, mat);
    box->SetPos(ChVector3d(-1, 0, 0.2));
    sys.Add(box);

    double W = box->GetMass() * g;
    for (int k = 0; k < 100; k++) {
        sys.DoStepDynamics(5e-3);
        EXPECT_LT(box->GetPosDt().Length(), 1e-5) << "step " << k;
        EXPECT_LT(box->GetAngVelParent().Length(), 1e-4) << "step " << k;
        EXPECT_NEAR(box->GetContactForce().z(), W, 1e-3 * W) << "step " << k;
    }
}
