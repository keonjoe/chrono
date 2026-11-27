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
// This project simulates the collision of a sphere against a fixed wall whose mass is
// smaller than the sphere's, and compares the measured coefficient of restitution against
// the input coefficient of restitution. A fixed body acts as one of infinite mass, so its
// mass must not enter the contact force.
//
// =============================================================================

#include "gtest/gtest.h"

#define SMC_MULTICORE
#include "../utest_SMC.h"

// Test system parameterized by SMC contact force model and input COR
class CorWallTest : public ::testing::TestWithParam<std::tuple<ChSystemSMC::ContactForceModel, float>> {
  protected:
    CorWallTest() {
        auto fmodel = std::get<0>(GetParam());

        auto mat = chrono_types::make_shared<ChContactMaterialSMC>();
        mat->SetYoungModulus(2.0e5f);
        mat->SetPoissonRatio(0.3f);
        mat->SetFriction(0.3f);
        mat->SetRestitution(std::get<1>(GetParam()));

        sys = new ChSystemMulticoreSMC();
        time_step = 3.0E-5;
        SetSimParameters(sys, ChVector3d(0, 0, 0), fmodel);

        sys->SetNumThreads(2);

        // Fixed wall with a mass of 0.1 (top face at y = 0)
        AddWall(sys, mat, ChVector3d(4, 1, 4), 0.1, ChVector3d(0, -0.5, 0), ChVector3d(0, 0, 0), true);

        // Sphere of mass 1 approaching the wall
        double srad = 0.5;
        ChVector3d init_v(0, -1, 0);
        body = AddSphere(sys, mat, srad, 1.0, ChVector3d(0, srad * 1.25, 0), init_v);
        vel_in = init_v.Length();
    }

    ~CorWallTest() { delete sys; }

    ChSystemMulticoreSMC* sys;
    std::shared_ptr<ChBody> body;
    double time_step;
    double vel_in;
};

TEST_P(CorWallTest, impact) {
    double t_end = 0.5;
    while (sys->GetChTime() < t_end) {
        sys->DoStepDynamics(time_step);
    }

    double cor_out = body->GetPosDt().y() / vel_in;
    double cor_in = std::get<1>(GetParam());
    std::cout << ForceModel_name(std::get<0>(GetParam())) << "  " << cor_in << "  " << cor_out << "\n";
    ASSERT_NEAR(cor_out, cor_in, 2e-3);
}

INSTANTIATE_TEST_SUITE_P(ChronoMulticore,
                         CorWallTest,
                         ::testing::Combine(::testing::Values(ChSystemSMC::ContactForceModel::Hooke, ChSystemSMC::ContactForceModel::Hertz), ::testing::Values(0.25, 0.5, 0.75)));
