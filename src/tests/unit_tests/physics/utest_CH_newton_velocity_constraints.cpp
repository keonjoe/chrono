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
// Test for constraints without position-level residual (ChShaftsGearbox, ChShaftsGear without phase-drift
// avoidance, ChShaftBodyRotation) with the Newton-iterated implicit integrators, when the velocities at the beginning
// of the step violate the constraint.
//
// Shaft A (inertia JA, speed w0) is coupled with ratio r to a second item of inertia JB (a shaft, or the rotation of a
// body about x), whose speed is reset to 0 before every step. The step must remove the mismatch with an impulse P:
//    JA*wA' = JA*wA - r*P,   JB*wB' = P,   wB' = r*wA'
// so that  wA' = wA * JA / (JA + r^2*JB),  P = JB*wB'  and, after n steps,  wA = w0 * (JA / (JA + r^2*JB))^n.
// The reaction reported by the integrators satisfies P = h*c*reaction, with c = 1 for implicit Euler and trapezoidal,
// c = gamma for Newmark and c = (1+alpha)*gamma for HHT (their end-of-step multiplier).
//
// =============================================================================

#include <cmath>
#include <functional>

#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChShaft.h"
#include "chrono/physics/ChShaftBodyConstraint.h"
#include "chrono/physics/ChShaftsGear.h"
#include "chrono/physics/ChShaftsGearbox.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/solver/ChDirectSolverLS.h"
#include "chrono/timestepper/ChTimestepperHHT.h"
#include "chrono/timestepper/ChTimestepperImplicit.h"

#include "gtest/gtest.h"

using namespace chrono;

enum class VelocityConstraint { GEARBOX, GEAR_NO_PHASE, SHAFT_BODY_ROTATION };

class NewtonVelocityConstraintTest : public ::testing::TestWithParam<std::tuple<ChTimestepper::Type, VelocityConstraint>> {};

TEST_P(NewtonVelocityConstraintTest, speed_reset) {
    ChTimestepper::Type type = std::get<0>(GetParam());
    VelocityConstraint constraint = std::get<1>(GetParam());

    double JA = 1;
    double JB = 0.01;
    double w0 = 100;
    double h = 1e-3;
    int num_steps = 100;

    ChSystemSMC sys;
    sys.SetGravitationalAcceleration(VNULL);
    sys.SetSolver(chrono_types::make_shared<ChSolverSparseQR>());
    switch (type) {
        case ChTimestepper::Type::EULER_IMPLICIT:
            sys.SetTimestepper(chrono_types::make_shared<ChTimestepperEulerImplicit>(&sys));
            break;
        case ChTimestepper::Type::TRAPEZOIDAL:
            sys.SetTimestepper(chrono_types::make_shared<ChTimestepperTrapezoidal>(&sys));
            break;
        case ChTimestepper::Type::NEWMARK:
            sys.SetTimestepper(chrono_types::make_shared<ChTimestepperNewmark>(&sys));
            break;
        case ChTimestepper::Type::HHT:
            sys.SetTimestepper(chrono_types::make_shared<ChTimestepperHHT>(&sys));
            break;
        default:
            sys.SetTimestepper(chrono_types::make_shared<ChTimestepperEulerImplicitLinearized>(&sys));
            break;
    }

    // impulse factor of the reaction
    double c = 1;
    if (auto hht = std::dynamic_pointer_cast<ChTimestepperHHT>(sys.GetTimestepper())) {
        double alpha = hht->GetAlpha();
        c = (1 + alpha) * (0.5 - alpha);
    } else if (auto newmark = std::dynamic_pointer_cast<ChTimestepperNewmark>(sys.GetTimestepper())) {
        c = newmark->GetGamma();
    }

    auto shaftA = chrono_types::make_shared<ChShaft>();
    shaftA->SetInertia(JA);
    shaftA->SetPosDt(w0);
    sys.AddShaft(shaftA);

    auto ground = chrono_types::make_shared<ChBody>();
    ground->SetFixed(true);
    sys.AddBody(ground);

    double r = 0.5;
    std::shared_ptr<ChShaft> shaftB;
    std::shared_ptr<ChBody> body;
    std::function<double()> reaction;
    switch (constraint) {
        case VelocityConstraint::GEARBOX: {
            shaftB = chrono_types::make_shared<ChShaft>();
            shaftB->SetInertia(JB);
            sys.AddShaft(shaftB);
            auto gearbox = chrono_types::make_shared<ChShaftsGearbox>();
            ChVector3d dir(1, 0, 0);
            gearbox->Initialize(shaftA, shaftB, ground, dir);
            gearbox->SetTransmissionRatio(r);
            sys.Add(gearbox);
            reaction = [gearbox]() { return gearbox->GetReaction2(); };
            break;
        }
        case VelocityConstraint::GEAR_NO_PHASE: {
            shaftB = chrono_types::make_shared<ChShaft>();
            shaftB->SetInertia(JB);
            sys.AddShaft(shaftB);
            auto gear = chrono_types::make_shared<ChShaftsGear>();
            gear->Initialize(shaftA, shaftB);
            gear->SetTransmissionRatio(r);
            gear->AvoidPhaseDrift(false);
            sys.Add(gear);
            reaction = [gear]() { return gear->GetReaction2(); };
            break;
        }
        case VelocityConstraint::SHAFT_BODY_ROTATION: {
            r = 1;
            body = chrono_types::make_shared<ChBody>();
            body->SetMass(1);
            body->SetInertiaXX(ChVector3d(JB, 1, 1));
            sys.AddBody(body);
            auto link = chrono_types::make_shared<ChShaftBodyRotation>();
            ChVector3d dir(1, 0, 0);
            link->Initialize(shaftA, body, dir);
            sys.Add(link);
            reaction = [link]() { return -link->GetTorqueReactionOnShaft(); };
            break;
        }
    }

    auto speedB = [&]() { return shaftB ? shaftB->GetPosDt() : body->GetAngVelLocal().x(); };
    auto resetB = [&]() {
        if (shaftB)
            shaftB->SetPosDt(0);
        else
            body->SetAngVelLocal(VNULL);
    };

    double k = JA / (JA + r * r * JB);

    for (int i = 0; i < num_steps; i++) {
        resetB();
        sys.DoStepDynamics(h);

        if (i == 0) {
            double wA = w0 * k;
            double P = JB * r * wA;
            EXPECT_NEAR(shaftA->GetPosDt(), wA, 1e-9 * w0);
            EXPECT_NEAR(speedB(), r * wA, 1e-9 * w0);
            EXPECT_NEAR(h * c * std::abs(reaction()), P, 1e-6 * P);
        }
    }

    EXPECT_NEAR(shaftA->GetPosDt(), w0 * std::pow(k, num_steps), 1e-8 * w0);
    EXPECT_NEAR(speedB(), r * shaftA->GetPosDt(), 1e-9 * w0);
}

static const std::tuple<ChTimestepper::Type, VelocityConstraint> velocity_constraint_cases[] = {
    {ChTimestepper::Type::EULER_IMPLICIT, VelocityConstraint::GEARBOX},
    {ChTimestepper::Type::EULER_IMPLICIT, VelocityConstraint::GEAR_NO_PHASE},
    {ChTimestepper::Type::EULER_IMPLICIT, VelocityConstraint::SHAFT_BODY_ROTATION},
    {ChTimestepper::Type::TRAPEZOIDAL, VelocityConstraint::GEARBOX},
    {ChTimestepper::Type::NEWMARK, VelocityConstraint::GEARBOX},
    {ChTimestepper::Type::NEWMARK, VelocityConstraint::SHAFT_BODY_ROTATION},
    {ChTimestepper::Type::HHT, VelocityConstraint::GEARBOX},
    {ChTimestepper::Type::HHT, VelocityConstraint::GEAR_NO_PHASE},
    {ChTimestepper::Type::HHT, VelocityConstraint::SHAFT_BODY_ROTATION},
    {ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED, VelocityConstraint::GEARBOX},
};

INSTANTIATE_TEST_SUITE_P(Physics, NewtonVelocityConstraintTest, ::testing::ValuesIn(velocity_constraint_cases));
