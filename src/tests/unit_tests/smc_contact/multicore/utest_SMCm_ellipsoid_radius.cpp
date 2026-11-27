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
// This project simulates an ellipsoid resting under its weight on a fixed wall and
// compares the penetration against the Hertz solution. The Hertz contact stiffness
// depends on the effective radius of curvature of the contact, which for an ellipsoid
// with semi-axes (a, b, c) resting on its c axis on a flat face is a * b / c. The
// effective radius must therefore follow the size of the shape.
//
// =============================================================================

#include "gtest/gtest.h"

#define SMC_MULTICORE
#include "../utest_SMC.h"

// Test system parameterized by the size (scale factor) of the ellipsoid
class EllipsoidRadiusTest : public ::testing::TestWithParam<double> {
  protected:
    EllipsoidRadiusTest() {
        double scale = GetParam();

        float young = 2.0e6f;
        float poisson = 0.3f;

        auto mat = chrono_types::make_shared<ChContactMaterialSMC>();
        mat->SetYoungModulus(young);
        mat->SetPoissonRatio(poisson);
        mat->SetFriction(0.3f);
        mat->SetRestitution(0.01f);

        sys = new ChSystemMulticoreSMC();
        time_step = 1.0e-4;
        double gravity = 9.81;
        SetSimParameters(sys, ChVector3d(0, -gravity, 0), ChSystemSMC::ContactForceModel::Hertz);

        sys->SetNumThreads(2);

        // Fixed wall (top face at y = 0)
        AddWall(sys, mat, scale * ChVector3d(4, 1, 4), 1.0, scale * ChVector3d(0, -0.5, 0), ChVector3d(0, 0, 0), true);

        // Ellipsoid touching the wall, with its shortest semi-axis vertical
        ChVector3d semiaxes = scale * ChVector3d(0.04, 0.02, 0.03);
        double mass = 1000 * CH_4_3 * CH_PI * semiaxes.x() * semiaxes.y() * semiaxes.z();
        ChVector3d inertia(0.2 * mass * (semiaxes.y() * semiaxes.y() + semiaxes.z() * semiaxes.z()), 0.2 * mass * (semiaxes.x() * semiaxes.x() + semiaxes.z() * semiaxes.z()),
                           0.2 * mass * (semiaxes.x() * semiaxes.x() + semiaxes.y() * semiaxes.y()));

        body = chrono_types::make_shared<ChBody>();
        body->SetMass(mass);
        body->SetInertiaXX(inertia);
        body->SetPos(ChVector3d(0, semiaxes.y(), 0));
        body->EnableCollision(true);
        utils::AddEllipsoidGeometry(body.get(), mat, 2.0 * semiaxes);
        sys->AddBody(body);

        height = semiaxes.y();

        // Hertz penetration, F = 4/3 * E_eff * sqrt(R_eff) * delta^(3/2)
        double E_eff = young / (2 * (1 - poisson * poisson));
        double R_eff = semiaxes.x() * semiaxes.z() / semiaxes.y();
        delta = std::pow(3 * mass * gravity / (4 * E_eff * std::sqrt(R_eff)), 2.0 / 3);
    }

    ~EllipsoidRadiusTest() { delete sys; }

    ChSystemMulticoreSMC* sys;
    std::shared_ptr<ChBody> body;
    double time_step;
    double height;
    double delta;
};

TEST_P(EllipsoidRadiusTest, penetration) {
    double t_end = 1.5;
    while (sys->GetChTime() < t_end) {
        sys->DoStepDynamics(time_step);
    }

    double penetration = height - body->GetPos().y();
    std::cout << "scale " << GetParam() << "  penetration " << penetration << "  Hertz " << delta << "\n";
    ASSERT_LT(std::abs(body->GetPosDt().y()), 1e-4 * GetParam());
    ASSERT_NEAR(penetration, delta, 0.02 * delta);
}

INSTANTIATE_TEST_SUITE_P(ChronoMulticore, EllipsoidRadiusTest, ::testing::Values(1.0, 10.0));
