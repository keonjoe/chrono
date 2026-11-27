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
// Test for the effective radius of curvature of contacts reported by the Bullet
// collision system, used by the Hertzian SMC contact force models.
//
// - indentation: a sphere of radius R rests on a box. With the Hertzian force
//   F = 4/3 E* sqrt(R) delta^(3/2), the static indentation under its weight F = m g is
//   delta = (9 F^2 / (16 E*^2 R))^(1/3).
// - radius: effective radius 1 / (1 / R1 + 1 / R2) for a sphere on a sphere, an ellipsoid
//   on a plane (a * b / c at the tip of the c semi-axis) and a cylinder lying on its side.
//
// =============================================================================

#include <cmath>

#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/collision/bullet/ChCollisionSystemBullet.h"

#include "gtest/gtest.h"

using namespace chrono;

class EffRadiusReporter : public ChContactContainer::ReportContactCallback {
  public:
    virtual bool OnReportContact(const ChVector3d& pA,
                                 const ChVector3d& pB,
                                 const ChMatrix33<>& plane_coord,
                                 double distance,
                                 double eff_radius,
                                 const ChVector3d& react_forces,
                                 const ChVector3d& react_torques,
                                 ChContactable* contactobjA,
                                 ChContactable* contactobjB,
                                 int constraint_offset) override {
        radius = eff_radius;
        num_contacts++;
        return true;
    }
    double radius = 0;
    int num_contacts = 0;
};

static std::shared_ptr<ChContactMaterialSMC> Material() {
    auto mat = chrono_types::make_shared<ChContactMaterialSMC>();
    mat->SetYoungModulus(2e7f);
    mat->SetPoissonRatio(0.3f);
    mat->SetRestitution(0.2f);
    mat->SetFriction(0.5f);
    return mat;
}

// Create an SMC system with a fixed box (top face at z = 0) and the given body
static void CreateSystem(ChSystemSMC& sys, std::shared_ptr<ChBody> body) {
    sys.SetCollisionSystem(chrono_types::make_shared<ChCollisionSystemBullet>());
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, -9.81));
    sys.SetContactForceModel(ChSystemSMC::Hertz);

    auto ground = chrono_types::make_shared<ChBodyEasyBox>(4, 4, 0.2, 1000, false, true, Material());
    ground->SetPos(ChVector3d(0, 0, -0.1));
    ground->SetFixed(true);
    sys.Add(ground);
    sys.Add(body);
}

// Effective radius of the contact(s) after a few steps
static double EffectiveRadius(ChSystemSMC& sys) {
    for (int i = 0; i < 20; i++)
        sys.DoStepDynamics(1e-4);
    auto reporter = chrono_types::make_shared<EffRadiusReporter>();
    sys.GetContactContainer()->ReportAllContacts(reporter);
    EXPECT_GT(reporter->num_contacts, 0);
    return reporter->radius;
}

TEST(BulletEffectiveRadius, indentation) {
    for (double R : {0.03, 0.3}) {
        auto ball = chrono_types::make_shared<ChBodyEasySphere>(R, 1000, false, true, Material());
        ball->SetPos(ChVector3d(0, 0, R));
        ChSystemSMC sys;
        CreateSystem(sys, ball);

        while (sys.GetChTime() < 1.0)
            sys.DoStepDynamics(1e-4);

        double E_eff = 2e7 / (2 * (1 - 0.3 * 0.3));
        double F = ball->GetMass() * 9.81;
        double delta_ref = std::cbrt(9 * F * F / (16 * E_eff * E_eff * R));
        double delta = R - ball->GetPos().z();
        std::cout << "R = " << R << "  indentation = " << delta << "  Hertz = " << delta_ref << std::endl;
        EXPECT_NEAR(delta, delta_ref, 1e-3 * delta_ref);
    }
}

TEST(BulletEffectiveRadius, radius) {
    {
        // Sphere on a fixed sphere of the same radius
        ChSystemSMC sys;
        auto ball = chrono_types::make_shared<ChBodyEasySphere>(0.1, 1000, false, true, Material());
        ball->SetPos(ChVector3d(0, 0, 0.2 + 0.1));
        CreateSystem(sys, ball);
        auto base = chrono_types::make_shared<ChBodyEasySphere>(0.1, 1000, false, true, Material());
        base->SetPos(ChVector3d(0, 0, 0.1));
        base->SetFixed(true);
        sys.Add(base);
        EXPECT_NEAR(EffectiveRadius(sys), 0.05, 1e-6);
    }
    {
        // Ellipsoid with semi-axes (0.2, 0.1, 0.05) resting on the tip of its z semi-axis
        ChSystemSMC sys;
        auto ell = chrono_types::make_shared<ChBodyEasyEllipsoid>(ChVector3d(0.4, 0.2, 0.1), 1000, false, true, Material());
        ell->SetPos(ChVector3d(0, 0, 0.05));
        CreateSystem(sys, ell);
        EXPECT_NEAR(EffectiveRadius(sys), 0.2 * 0.1 / 0.05, 1e-3);
    }
    {
        // Cylinder of radius 0.05 lying on its side
        ChSystemSMC sys;
        auto cyl = chrono_types::make_shared<ChBodyEasyCylinder>(ChAxis::Z, 0.05, 0.3, 1000, false, true, Material());
        cyl->SetPos(ChVector3d(0, 0, 0.05));
        cyl->SetRot(QuatFromAngleX(CH_PI_2));
        CreateSystem(sys, cyl);
        EXPECT_NEAR(EffectiveRadius(sys), 0.05, 1e-6);
    }
}
