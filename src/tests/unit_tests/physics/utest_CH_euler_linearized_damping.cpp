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
// Test for the treatment of damping forces in the linearized implicit Euler timestepper.
//
// A free body of mass m is attached to ground by a linear damper (no spring, coefficient c)
// and given an initial velocity v0 along the damper axis. The analytical solution is
//    v(t) = v0 * exp(-c/m * t)
//    x(t) = x0 + v0 * m/c * (1 - exp(-c/m * t))
// With a damping Jacobian (stiff damper), the linearized Euler step is the implicit step
//    v_(n+1) = v_n * m / (m + c*dt)
// Without a damping Jacobian, the damping force is treated explicitly:
//    v_(n+1) = v_n * (1 - c*dt/m)
// In both cases x_(n+1) = x_n + dt * v_(n+1), and the solution converges to the analytical one as dt -> 0.
//
// The same damper is also integrated as a user-defined ChIntegrableIIorder, which relies on the default
// (finite difference) implementations of dF/dv.
//
// =============================================================================

#include "chrono/physics/ChLinkTSDA.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/solver/ChDirectSolverLS.h"
#include "chrono/solver/ChIterativeSolverLS.h"
#include "chrono/solver/ChSolverPSOR.h"
#include "chrono/timestepper/ChTimestepperImplicit.h"

#include "gtest/gtest.h"

using namespace chrono;

enum class DamperSolver { SPARSE_QR, MINRES, PSOR };

class LinearizedEulerDampingTest : public ::testing::TestWithParam<std::tuple<DamperSolver, bool>> {};

TEST_P(LinearizedEulerDampingTest, linear_damper) {
    DamperSolver solver_type = std::get<0>(GetParam());
    bool stiff = std::get<1>(GetParam());

    double m = 1;   // body mass
    double c = 2;   // damping coefficient
    double v0 = 1;  // initial velocity
    double x0 = 1;  // initial position
    double dt = 1e-3;
    double t_end = 1;

    ChSystemNSC sys;
    sys.SetGravitationalAcceleration(VNULL);
    sys.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);
    switch (solver_type) {
        case DamperSolver::SPARSE_QR:
            sys.SetSolver(chrono_types::make_shared<ChSolverSparseQR>());
            break;
        case DamperSolver::MINRES: {
            auto solver = chrono_types::make_shared<ChSolverMINRES>();
            solver->SetMaxIterations(100);
            solver->SetTolerance(1e-14);
            sys.SetSolver(solver);
            break;
        }
        case DamperSolver::PSOR:
            sys.SetSolver(chrono_types::make_shared<ChSolverPSOR>());
            break;
    }

    auto ground = chrono_types::make_shared<ChBody>();
    ground->SetFixed(true);
    sys.Add(ground);

    auto body = chrono_types::make_shared<ChBody>();
    body->SetMass(m);
    body->SetInertiaXX(ChVector3d(1, 1, 1));
    body->SetPos(ChVector3d(x0, 0, 0));
    body->SetPosDt(ChVector3d(v0, 0, 0));
    sys.Add(body);

    auto damper = chrono_types::make_shared<ChLinkTSDA>();
    damper->Initialize(ground, body, true, VNULL, VNULL);
    damper->SetSpringCoefficient(0);
    damper->SetDampingCoefficient(c);
    damper->IsStiff(stiff);
    sys.AddLink(damper);

    // velocity amplification factor of one step (implicit with a damping Jacobian, explicit otherwise)
    double factor = stiff ? m / (m + c * dt) : 1 - c * dt / m;

    // discrete solution: v_(n+1) = factor * v_n,  x_(n+1) = x_n + dt * v_(n+1)
    double v_discrete = v0;
    double x_discrete = x0;
    while (sys.GetChTime() < t_end - dt / 2) {
        sys.DoStepDynamics(dt);
        v_discrete *= factor;
        x_discrete += dt * v_discrete;
    }

    double v = body->GetPosDt().x();
    double x = body->GetPos().x();

    double v_exact = v0 * std::exp(-c / m * t_end);
    double x_exact = x0 + v0 * m / c * (1 - std::exp(-c / m * t_end));

    EXPECT_NEAR(v, v_discrete, 1e-8);
    EXPECT_NEAR(x, x_discrete, 1e-8);

    // first order convergence to the analytical solution
    EXPECT_NEAR(v, v_exact, 2 * c * dt / m * v0);
    EXPECT_NEAR(x, x_exact, 2 * c * dt / m * v0 * m / c);
    EXPECT_NEAR(body->GetPosDt().y(), 0, 1e-12);
    EXPECT_NEAR(body->GetPosDt().z(), 0, 1e-12);
}

// VI solvers (PSOR) do not support stiffness or damping matrices, so only the non-stiff damper is tested with them
INSTANTIATE_TEST_SUITE_P(Physics,
                         LinearizedEulerDampingTest,
                         ::testing::Values(std::make_tuple(DamperSolver::SPARSE_QR, true),
                                           std::make_tuple(DamperSolver::SPARSE_QR, false),
                                           std::make_tuple(DamperSolver::MINRES, true),
                                           std::make_tuple(DamperSolver::MINRES, false),
                                           std::make_tuple(DamperSolver::PSOR, false)));

// One-dimensional damper  m*a = -c*v, defined only through LoadResidual_F and LoadResidual_Mv
class DamperIntegrable : public ChIntegrableIIorder {
  public:
    DamperIntegrable(double m, double c, double x0, double v0) : m_m(m), m_c(c), m_x(x0), m_v(v0), m_T(0) {}

    virtual unsigned int GetNumCoordsPosLevel() override { return 1; }
    virtual unsigned int GetNumCoordsVelLevel() override { return 1; }
    virtual unsigned int GetNumCoordsAccLevel() override { return 1; }

    virtual void StateGather(ChState& x, ChStateDelta& v, double& T) override {
        x(0) = m_x;
        v(0) = m_v;
        T = m_T;
    }

    virtual void StateScatter(const ChState& x, const ChStateDelta& v, const double T, UpdateFlags update_flags) override {
        m_x = x(0);
        m_v = v(0);
        m_T = T;
    }

    virtual void LoadResidual_F(ChVectorDynamic<>& R, const double c) override { R(0) += c * (-m_c * m_v); }

    virtual void LoadResidual_Mv(ChVectorDynamic<>& R, const ChVectorDynamic<>& w, const double c) override { R(0) += c * m_m * w(0); }

  private:
    double m_m;
    double m_c;
    double m_x;
    double m_v;
    double m_T;
};

TEST(LinearizedEulerDamping, integrable) {
    double m = 1;   // mass
    double c = 2;   // damping coefficient
    double v0 = 1;  // initial velocity
    double x0 = 1;  // initial position
    double dt = 1e-3;
    double t_end = 1;

    DamperIntegrable integrable(m, c, x0, v0);
    ChTimestepperEulerImplicitLinearized stepper(&integrable);

    double v_discrete = v0;
    while (stepper.GetTime() < t_end - dt / 2) {
        stepper.Advance(dt);
        v_discrete *= m / (m + c * dt);
    }

    double v_exact = v0 * std::exp(-c / m * t_end);

    // the finite difference Jacobian is exact up to round-off for a linear damper
    EXPECT_NEAR(stepper.GetStateVel()(0), v_discrete, 1e-6);
    EXPECT_NEAR(stepper.GetStateVel()(0), v_exact, 2 * c * dt / m * v0);
}
