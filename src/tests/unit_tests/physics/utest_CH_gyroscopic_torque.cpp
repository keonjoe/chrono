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
// Test for the gyroscopic torque with the integrators that evaluate the forces once per step.
//
// 1. A torque-free symmetric body (inertia I1 = I2, I3) spins with w3 about its symmetry axis plus a transverse
//    component w_t. In the body frame, w3 is constant and the transverse component rotates at the nutation rate
//    W = (I1 - I3) / I1 * w3, so the kinetic energy and |L| are constant. With an explicit gyroscopic torque the
//    nutation grows by sqrt(1 + (h*W)^2) per step; with the linearly implicit (Kahan) torque, E and |L| are conserved
//    for any step and the transverse velocity turns by 2*atan(h*W/2) per step.
// 2. A torque-free triaxial body spinning close to its intermediate axis flips with the period 4*K(k)/lambda
//    (Landau & Lifshitz, Mechanics, sec. 37).
// 3. A heavy symmetric top on a spherical joint (pivot at distance l from the center of mass), started in steady
//    precession at tilt theta: the precession rate p satisfies m*g*l = p*(I3*w3 - I1p*p*cos(theta)), with I1p the
//    transverse inertia about the pivot, and the tilt stays constant. This checks that the gyroscopic torque balanced
//    by the joint reaction is not biased.
//
// =============================================================================

#include <algorithm>
#include <cmath>
#include <vector>

#include "chrono/physics/ChBodyAuxRef.h"
#include "chrono/physics/ChLinkLock.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/solver/ChDirectSolverLS.h"
#include "chrono/solver/ChSolverPSOR.h"
#include "chrono/timestepper/ChTimestepperExplicit.h"
#include "chrono/timestepper/ChTimestepperImplicit.h"

#include "gtest/gtest.h"

using namespace chrono;

enum class GyroSolver { SPARSE_QR, PSOR };

static void SetupSystem(ChSystemNSC& sys, ChTimestepper::Type stepper, GyroSolver solver) {
    if (stepper == ChTimestepper::Type::EULER_SEMI_IMPLICIT)
        sys.SetTimestepper(chrono_types::make_shared<ChTimestepperEulerSemiImplicit>(&sys));
    else
        sys.SetTimestepperType(stepper);
    if (solver == GyroSolver::SPARSE_QR) {
        sys.SetSolver(chrono_types::make_shared<ChSolverSparseQR>());
    } else {
        auto psor = chrono_types::make_shared<ChSolverPSOR>();
        psor->SetMaxIterations(100);
        sys.SetSolver(psor);
    }
}

class GyroscopicTorqueTest : public ::testing::TestWithParam<std::tuple<ChTimestepper::Type, GyroSolver, bool>> {};

// Torque-free symmetric body
TEST_P(GyroscopicTorqueTest, free_symmetric_body) {
    ChSystemNSC sys;
    sys.SetGravitationalAcceleration(VNULL);
    SetupSystem(sys, std::get<0>(GetParam()), std::get<1>(GetParam()));

    double I1 = 2, I3 = 1;
    double w3 = 10, wt = 2;
    double dt = 1e-2;
    double t_end = 10;

    std::shared_ptr<ChBody> body;
    if (std::get<2>(GetParam())) {
        auto body_auxref = chrono_types::make_shared<ChBodyAuxRef>();
        body_auxref->SetFrameCOMToRef(ChFrame<>(ChVector3d(0.1, -0.2, 0.3), QuatFromAngleX(0.5)));
        body = body_auxref;
    } else {
        body = chrono_types::make_shared<ChBody>();
    }
    body->SetMass(1);
    body->SetInertiaXX(ChVector3d(I1, I1, I3));
    sys.Add(body);
    body->SetAngVelLocal(ChVector3d(wt, 0, w3));

    ChMatrix33<> J(ChVector3d(I1, I1, I3));
    ChVector3d L0 = body->GetRotMat() * (J * body->GetAngVelLocal());
    double E0 = 0.5 * Vdot(body->GetAngVelLocal(), J * body->GetAngVelLocal());

    // nutation phase: angle of the transverse velocity in the body frame
    double W = (I1 - I3) / I1 * w3;
    double phase = 0;
    double phase_discrete = 0;
    ChVector3d w_old = body->GetAngVelLocal();
    while (sys.GetChTime() < t_end - dt / 2) {
        sys.DoStepDynamics(dt);
        ChVector3d w = body->GetAngVelLocal();
        phase += std::atan2(w_old.x() * w.y() - w_old.y() * w.x(), w_old.x() * w.x() + w_old.y() * w.y());
        phase_discrete += 2 * std::atan(dt * W / 2);
        w_old = w;
    }

    ChVector3d w = body->GetAngVelLocal();
    ChVector3d L = body->GetRotMat() * (J * w);
    double E = 0.5 * Vdot(w, J * w);

    // conservation of energy and of the magnitude of the angular momentum (explicit torque: E/E0 = 1.83)
    EXPECT_NEAR(E / E0, 1, 1e-9);
    EXPECT_NEAR(L.Length() / L0.Length(), 1, 1e-9);
    EXPECT_NEAR(w.z(), w3, 1e-9);
    EXPECT_NEAR(std::sqrt(w.x() * w.x() + w.y() * w.y()), wt, 1e-9);

    // nutation: the body-frame rotation of the transverse velocity is -W per unit time (Cayley angle per step)
    EXPECT_NEAR(phase, -phase_discrete, 1e-6);
    EXPECT_NEAR(phase, -W * t_end, 0.05);

    // direction of the angular momentum in the absolute frame (bounded error of order dt)
    EXPECT_LT(std::acos(std::min(1.0, Vdot(L, L0) / (L.Length() * L0.Length()))), 0.05);
}

// integrator, solver, body with a reference frame other than the center of mass
static const std::tuple<ChTimestepper::Type, GyroSolver, bool> gyro_cases[] = {
    {ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED, GyroSolver::SPARSE_QR, false}, {ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED, GyroSolver::PSOR, false},
    {ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED, GyroSolver::SPARSE_QR, true},  {ChTimestepper::Type::EULER_IMPLICIT_PROJECTED, GyroSolver::SPARSE_QR, false},
    {ChTimestepper::Type::TRAPEZOIDAL_LINEARIZED, GyroSolver::SPARSE_QR, false},    {ChTimestepper::Type::EULER_SEMI_IMPLICIT, GyroSolver::SPARSE_QR, false},
};

INSTANTIATE_TEST_SUITE_P(Physics, GyroscopicTorqueTest, ::testing::ValuesIn(gyro_cases));

// Complete elliptic integral of the first kind, by the arithmetic-geometric mean
static double EllipticK(double k) {
    double a = 1, b = std::sqrt(1 - k * k);
    for (int i = 0; i < 30; i++) {
        double a_next = 0.5 * (a + b);
        b = std::sqrt(a * b);
        a = a_next;
    }
    return CH_PI / (2 * a);
}

// Torque-free triaxial body close to its intermediate axis
TEST(GyroscopicTorque, intermediate_axis_flip) {
    ChSystemNSC sys;
    sys.SetGravitationalAcceleration(VNULL);
    SetupSystem(sys, ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED, GyroSolver::PSOR);

    double I1 = 1, I2 = 2, I3 = 3;
    ChVector3d w0(0.5, 10, 0.5);
    double dt = 2e-3;
    double t_end = 15;

    auto body = chrono_types::make_shared<ChBody>();
    body->SetMass(1);
    body->SetInertiaXX(ChVector3d(I1, I2, I3));
    sys.Add(body);
    body->SetAngVelLocal(w0);

    // period of the flip (L^2 > 2*E*I2)
    double E2 = I1 * w0.x() * w0.x() + I2 * w0.y() * w0.y() + I3 * w0.z() * w0.z();
    double L2 = I1 * I1 * w0.x() * w0.x() + I2 * I2 * w0.y() * w0.y() + I3 * I3 * w0.z() * w0.z();
    ASSERT_GT(L2, E2 * I2);
    double lambda = std::sqrt((I3 - I2) * (L2 - E2 * I1) / (I1 * I2 * I3));
    double k2 = (I2 - I1) * (E2 * I3 - L2) / ((I3 - I2) * (L2 - E2 * I1));
    double period = 4 * EllipticK(std::sqrt(k2)) / lambda;

    // times at which w_y changes sign from negative to positive
    std::vector<double> crossings;
    double wy_old = w0.y();
    double E_err = 0;
    while (sys.GetChTime() < t_end - dt / 2) {
        sys.DoStepDynamics(dt);
        ChVector3d w = body->GetAngVelLocal();
        if (wy_old < 0 && w.y() >= 0)
            crossings.push_back(sys.GetChTime() - dt * w.y() / (w.y() - wy_old));
        wy_old = w.y();
        double E = I1 * w.x() * w.x() + I2 * w.y() * w.y() + I3 * w.z() * w.z();
        E_err = std::max(E_err, std::abs(E / E2 - 1));
    }

    // explicit torque: the energy error reaches 31% and the period is off by 1e-3
    ASSERT_GE(crossings.size(), 3u);
    double period_num = (crossings.back() - crossings.front()) / (crossings.size() - 1);
    EXPECT_NEAR(period_num / period, 1, 1e-4);
    EXPECT_LT(E_err, 1e-3);
}

// Heavy symmetric top on a spherical joint, in steady precession
TEST(GyroscopicTorque, heavy_top_steady_precession) {
    for (auto solver : {GyroSolver::SPARSE_QR, GyroSolver::PSOR}) {
        ChSystemNSC sys;
        double g = 9.81;
        sys.SetGravitationalAcceleration(ChVector3d(0, 0, -g));
        SetupSystem(sys, ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED, solver);

        double m = 1;
        double I1 = 2e-3, I3 = 1e-3;  // central inertia
        double l = 0.05;              // distance from the pivot to the center of mass
        double theta = 10 * CH_DEG_TO_RAD;
        double w3 = 300;
        double dt = 1e-4;
        double t_end = 1;

        // steady precession rate (slow solution)
        double I1p = I1 + m * l * l;
        double c = std::cos(theta);
        double p = (I3 * w3 - std::sqrt(I3 * I3 * w3 * w3 - 4 * I1p * c * m * g * l)) / (2 * I1p * c);

        auto ground = chrono_types::make_shared<ChBody>();
        ground->SetFixed(true);
        sys.Add(ground);

        auto top = chrono_types::make_shared<ChBody>();
        top->SetMass(m);
        top->SetInertiaXX(ChVector3d(I1, I1, I3));
        ChQuaternion<> q0 = QuatFromAngleX(theta);
        top->SetRot(q0);
        top->SetPos(top->GetRotMat() * ChVector3d(0, 0, l));
        // w = p * (vertical axis) + psi_dot * (symmetry axis), with w3 the symmetry axis component
        ChVector3d vertical_local = top->GetRotMat().transpose() * ChVector3d(0, 0, 1);
        ChVector3d w_local = p * vertical_local + ChVector3d(0, 0, w3 - p * vertical_local.z());
        top->SetAngVelLocal(w_local);
        top->SetPosDt(Vcross(top->GetRotMat() * w_local, top->GetPos()));
        sys.Add(top);

        auto joint = chrono_types::make_shared<ChLinkLockSpherical>();
        joint->Initialize(top, ground, ChFrame<>(VNULL));
        sys.AddLink(joint);

        double azimuth = 0;
        double azimuth_old = std::atan2(top->GetRotMat().GetAxisZ().y(), top->GetRotMat().GetAxisZ().x());
        double tilt_err = 0;
        while (sys.GetChTime() < t_end - dt / 2) {
            sys.DoStepDynamics(dt);
            ChVector3d axis = top->GetRotMat().GetAxisZ();
            double a = std::atan2(axis.y(), axis.x());
            azimuth += std::remainder(a - azimuth_old, CH_2PI);
            azimuth_old = a;
            tilt_err = std::max(tilt_err, std::abs(std::acos(axis.z()) - theta));
        }

        // tilt error with the explicit torque: 0.26 deg; with a predictor that ignores the joint reaction: 0.13 deg
        EXPECT_LT(tilt_err, 0.02 * CH_DEG_TO_RAD);
        EXPECT_NEAR(azimuth / t_end, p, 2e-3 * p);
    }
}
