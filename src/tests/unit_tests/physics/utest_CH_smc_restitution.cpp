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
// Test for the coefficient of restitution of the default SMC contact force models.
//
// A sphere approaches a fixed box with speed v (no gravity). The normal damping of the Hooke
// and Hertz models is derived from the material coefficient of restitution e, and the
// contact force is never attractive, so the objects separate as soon as it vanishes. The
// rebound speed must be e * v, for any impact speed (both models have a velocity-independent
// coefficient of restitution with material-based coefficients).
//
// - impact: fixed box of large mass, several values of e and of the impact speed.
// - light_ground: the fixed box has a mass smaller than the sphere's. A fixed body acts as
//   one of infinite mass, so the result must be the same.
//
// =============================================================================

#include <cmath>

#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemSMC.h"

#include "gtest/gtest.h"

using namespace chrono;

// Returns the ratio of rebound to impact speed
static double Rebound(ChSystemSMC::ContactForceModel model, double e, double v, double ground_mass) {
    double radius = 0.1;
    double E = 2e7;
    double nu = 0.3;

    ChSystemSMC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, 0));
    sys.SetContactForceModel(model);
    sys.SetCharacteristicImpactVelocity(v);

    auto mat = chrono_types::make_shared<ChContactMaterialSMC>();
    mat->SetYoungModulus((float)E);
    mat->SetPoissonRatio((float)nu);
    mat->SetFriction(0.5f);
    mat->SetRestitution((float)e);

    auto ground = chrono_types::make_shared<ChBodyEasyBox>(1, 1, 0.2, 1000, false, true, mat);
    ground->SetPos(ChVector3d(0, 0, -0.1));
    ground->SetFixed(true);
    if (ground_mass > 0)
        ground->SetMass(ground_mass);
    sys.Add(ground);

    auto ball = chrono_types::make_shared<ChBodyEasySphere>(radius, 1000, false, true, mat);
    ball->SetPos(ChVector3d(0, 0, radius + 1e-6));
    ball->SetPosDt(ChVector3d(0, 0, -v));
    sys.Add(ball);

    // Hertzian contact duration (Bullet uses the default effective radius of curvature) and step size
    double m = ball->GetMass();
    double E_eff = E / (2 * (1 - nu * nu));
    double R_eff = ChCollisionInfo::GetDefaultEffectiveCurvatureRadius();
    double tc = 2.87 * std::pow(m * m / (R_eff * E_eff * E_eff * v), 0.2);
    double step = tc / 200;

    bool touched = false;
    while (sys.GetChTime() < 20 * tc) {
        sys.DoStepDynamics(step);
        if (sys.GetNumContacts() > 0)
            touched = true;
        else if (touched)
            break;
    }
    EXPECT_TRUE(touched);

    return ball->GetPosDt().z() / v;
}

class SMCRestitutionTest : public ::testing::TestWithParam<ChSystemSMC::ContactForceModel> {};

TEST_P(SMCRestitutionTest, impact) {
    for (double e : {0.1, 0.2, 0.5, 0.8}) {
        for (double v : {0.5, 5.0}) {
            double cr = Rebound(GetParam(), e, v, 0);
            std::cout << "e = " << e << "  v = " << v << "  rebound/impact = " << cr << std::endl;
            EXPECT_NEAR(cr, e, 0.01 * e + 2e-3);
        }
    }
}

TEST_P(SMCRestitutionTest, light_ground) {
    double e = 0.5;
    double v = 2.0;
    double cr = Rebound(GetParam(), e, v, 0.1);
    std::cout << "e = " << e << "  v = " << v << "  rebound/impact = " << cr << std::endl;
    EXPECT_NEAR(cr, e, 0.01 * e + 2e-3);
}

INSTANTIATE_TEST_SUITE_P(ChSystemSMC, SMCRestitutionTest, ::testing::Values(ChSystemSMC::ContactForceModel::Hooke, ChSystemSMC::ContactForceModel::Hertz));
