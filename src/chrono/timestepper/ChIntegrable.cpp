// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Alessandro Tasora, Radu Serban
// =============================================================================

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <Eigen/LU>

#include "chrono/timestepper/ChIntegrable.h"

namespace chrono {

// -----------------------------------------------------------------------------
// Implementation of ChIntegrable methods
// -----------------------------------------------------------------------------

void ChIntegrable::StateIncrement(ChState& y_new, const ChState& y, const ChStateDelta& Dy) {
    //// RADU
    //// The next two lines don't make sense: we assert y and y_new are of the same size, then we resize y_new!!!

    assert(y_new.size() == y.size() && y.size() == Dy.size());

    //// RADU
    //// Fix this poor implementation!
    //// We cannot do the natural thing here (i.e, y_new = y + Dy) because this would lead to stack overflow!!!
    //// Indeed, see implementation of overloaded operator+()

    ////y_new = y + Dy;
    y_new.resize(y.size());
    for (int i = 0; i < y.size(); ++i) {
        y_new(i) = y(i) + Dy(i);
    }
}

// -----------------------------------------------------------------------------
// Implementation of ChIntegrableII methods
// -----------------------------------------------------------------------------

void ChIntegrableIIorder::StateSetup(ChState& x, ChStateDelta& v, ChStateDelta& a) {
    x.resize(GetNumCoordsPosLevel());
    v.resize(GetNumCoordsVelLevel());
    a.setZero(GetNumCoordsVelLevel(), a.GetIntegrable());  // defined value for items that do not store an acceleration
}

bool ChIntegrableIIorder::StateSolveA(ChStateDelta& Dvdt,        // result: computed a for a=dv/dt
                                      ChVectorDynamic<>& L,      // result: computed Lagrange multipliers
                                      const ChState& x,          // current state, x
                                      const ChStateDelta& v,     // current state, v
                                      const double T,            // current time T
                                      const double dt,           // timestep (if needed)
                                      bool force_state_scatter,  // if true, scatter x and v to the system
                                      UpdateFlags update_flags,  ///< if UpdateFlags::UPDATE_ALL, do a full update during scatter, otherwise switch off visual asset update, etc.
                                      ChLumpingParms* lumping    // use lumped masses to avoid inverting a mass matrix
) {
    if (force_state_scatter)
        StateScatter(x, v, T, update_flags);

    ChVectorDynamic<> R(GetNumCoordsVelLevel());

    if (!lumping) {
        ChVectorDynamic<> Qc(GetNumConstraints());
        const double Delta = 1e-6;

        R.setZero();
        Qc.setZero();

        LoadResidual_F(R, 1.0);
        LoadConstraint_C(Qc, -2.0 / (Delta * Delta), -1.0 / Delta);

        // numerical differentiation to get the Qc term in constraints
        ChStateDelta dx(v);
        dx *= Delta;
        ChState xdx(x.size(), this);

        StateIncrement(xdx, x, dx);
        StateScatter(xdx, v, T + Delta, update_flags);
        LoadConstraint_C(Qc, 1.0 / (Delta * Delta), 1.0 / Delta);

        StateIncrement(xdx, x, -dx);
        StateScatter(xdx, v, T - Delta, update_flags);
        LoadConstraint_C(Qc, 1.0 / (Delta * Delta), 0.0);

        StateScatter(x, v, T, update_flags);  // back to original state

        bool success = StateSolveCorrection(Dvdt, L, R, Qc, 1.0, 0, 0, x, v, T, false, UpdateFlags::UPDATE_ALL_NO_VISUAL, true, true);

        return success;
    } else {
        ChVectorDynamic<> Md(GetNumCoordsVelLevel());  // diagonal mass
        Md.setZero(GetNumCoordsVelLevel());

        // penalty terms from constraints
        R.setZero(GetNumCoordsVelLevel());
        L.setZero(GetNumConstraints());

        LoadResidual_F(R, 1.0);  // rhs, applied forces

        double err = 0;
        LoadLumpedMass_Md(Md, err, 1.0);

        if (GetNumConstraints()) {
            LoadConstraint_C(L, -lumping->Ck_penalty, -lumping->Ck_penalty);  // L  = -k*C     // to do: modulate  k  as constraint-dependent,
                                                        // k=lumping->Ck_penalty*Ck_i
            LoadResidual_CqL(R, L, 1.0);                // Fc =  Cq' * (-k*C)    = Cq' * L
            // to do: add c_penalty for speed proportional damping, as   Fc = - Cq' * (k*C + c*(dC/dt))
        }

        // compute acceleration using lumped diagonal mass. No need to invoke a linear solver.
        Dvdt.array() = R.array() / Md.array();  //  a = Md^-1 * F

        return (err == 0);
    }
}

void ChIntegrableIIorder::StateIncrementX(ChState& x_new, const ChState& x, const ChStateDelta& Dx) {
    //// RADU
    //// Fix this poor implementation!
    //// We cannot do the natural thing here (y_new = y + Dy) because this would lead to stack overflow!!!
    //// Indeed, see implementation of overloaded operator+()

    ////x_new = x + Dx;
    x_new.resize(x.size());
    for (int i = 0; i < x.size(); ++i) {
        x_new(i) = x(i) + Dx(i);
    }
}

void ChIntegrableIIorder::StateGather(ChState& y, double& T) {
    ChState mx(GetNumCoordsPosLevel(), y.GetIntegrable());
    ChStateDelta mv(GetNumCoordsVelLevel(), y.GetIntegrable());
    StateGather(mx, mv, T);
    y.segment(0, mx.size()) = mx;
    y.segment(GetNumCoordsPosLevel(), mv.size()) = mv;
}

void ChIntegrableIIorder::StateScatter(const ChState& y, const double T, UpdateFlags update_flags) {
    ChState mx(GetNumCoordsPosLevel(), y.GetIntegrable());
    ChStateDelta mv(GetNumCoordsVelLevel(), y.GetIntegrable());
    mx = y.segment(0, GetNumCoordsPosLevel());
    mv = y.segment(GetNumCoordsPosLevel(), GetNumCoordsVelLevel());
    StateScatter(mx, mv, T, update_flags);
}

void ChIntegrableIIorder::StateGatherDerivative(ChStateDelta& Dydt) {
    ChStateDelta mv(GetNumCoordsVelLevel(), Dydt.GetIntegrable());
    ChStateDelta ma(GetNumCoordsVelLevel(), Dydt.GetIntegrable());
    StateGatherAcceleration(ma);
    Dydt.segment(0, GetNumCoordsVelLevel()) = mv;
    Dydt.segment(GetNumCoordsVelLevel(), GetNumCoordsVelLevel()) = ma;
}

void ChIntegrableIIorder::StateScatterDerivative(const ChStateDelta& Dydt) {
    ChStateDelta ma(GetNumCoordsVelLevel(), Dydt.GetIntegrable());
    ma = Dydt.segment(GetNumCoordsVelLevel(), GetNumCoordsVelLevel());
    StateScatterAcceleration(ma);
}

void ChIntegrableIIorder::StateIncrement(ChState& y_new,         // resulting y_new = y + Dy
                                         const ChState& y,       // initial state y
                                         const ChStateDelta& Dy  // state increment Dy
) {
    if (y.size() == GetNumCoordsPosLevel()) {
        // Incrementing the x part only, user provided only x  in y={x, dx/dt}
        StateIncrementX(y_new, y, Dy);

        return;
    }

    if (y.size() == (GetNumCoordsPosLevel() + GetNumCoordsVelLevel())) {
        // Incrementing y in y={x, dx/dt}.
        // PERFORMANCE WARNING! temporary vectors allocated on heap.
        // This is only to support compatibility with 1st order integrators.
        ChState mx(GetNumCoordsPosLevel(), y.GetIntegrable());
        ChStateDelta mv(GetNumCoordsVelLevel(), y.GetIntegrable());
        mx = y.segment(0, GetNumCoordsPosLevel());
        mv = y.segment(GetNumCoordsPosLevel(), GetNumCoordsVelLevel());

        ChStateDelta mDx(GetNumCoordsVelLevel(), y.GetIntegrable());
        ChStateDelta mDv(GetNumCoordsVelLevel(), y.GetIntegrable());
        mDx = Dy.segment(0, GetNumCoordsVelLevel());
        mDv = Dy.segment(GetNumCoordsVelLevel(), GetNumCoordsVelLevel());

        ChState mx_new(GetNumCoordsPosLevel(), y.GetIntegrable());
        ChStateDelta mv_new(GetNumCoordsVelLevel(), y.GetIntegrable());
        StateIncrementX(mx_new, mx, mDx);  // increment positions
        mv_new = mv + mDv;                 // increment speeds

        y_new.segment(0, GetNumCoordsPosLevel()) = mx_new;
        y_new.segment(GetNumCoordsPosLevel(), GetNumCoordsVelLevel()) = mv_new;

        return;
    }

    throw std::runtime_error("StateIncrement() called with a wrong number of elements");
}

bool ChIntegrableIIorder::StateSolve(ChStateDelta& dydt,        // result: computed dydt
                                     ChVectorDynamic<>& L,      // result: computed Lagrange multipliers
                                     const ChState& y,          // current state y
                                     const double T,            // current time T
                                     const double dt,           // timestep (if needed, e.g. in NSC)
                                     bool force_state_scatter,  // if false, y and T are not scattered to the system
                                     UpdateFlags update_flags,   // if UpdateFlags::UPDATE_ALL, do a full update during scatter, otherwise switch off visual asset update, etc.
                                     ChLumpingParms* lumping    // if not null, uses lumped masses to avoid inverting a
                                                                // mass matrix, and uses penalty for constraints.
) {
    ChState mx(GetNumCoordsPosLevel(), y.GetIntegrable());
    ChStateDelta mv(GetNumCoordsVelLevel(), y.GetIntegrable());
    mx = y.segment(0, GetNumCoordsPosLevel());
    mv = y.segment(GetNumCoordsPosLevel(), GetNumCoordsVelLevel());

    // Solve with custom II order solver
    ChStateDelta ma(GetNumCoordsVelLevel(), y.GetIntegrable());
    if (!StateSolveA(ma, L, mx, mv, T, dt, force_state_scatter, update_flags, lumping)) {
        return false;
    }

    dydt.segment(0, GetNumCoordsVelLevel()) = mv;
    dydt.segment(GetNumCoordsVelLevel(), GetNumCoordsVelLevel()) = ma;

    return true;
}

// -----------------------------------------------------------------------------
// Default implementations of the ChIntegrable implicit integration interface
// -----------------------------------------------------------------------------

namespace {

// Finite difference step for a coordinate of given magnitude.
inline double FDStep(double val) {
    return 1e-7 * std::max(1.0, std::abs(val));
}

// Solve the (dense) Newton system
//   | G    CqT | |  D | = |  R |
//   | CqT' 0   | | -DL|   | -Qc|
// (only the first block row if there are no constraint columns). Returns false if the system matrix is singular.
bool SolveDenseKKT(const ChMatrixDynamic<>& G,
                   const ChMatrixDynamic<>& CqT,
                   const ChVectorDynamic<>& R,
                   const ChVectorDynamic<>& Qc,
                   ChVectorDynamic<>& D,
                   ChVectorDynamic<>& DL,
                   int num_constraints) {
    Eigen::Index n = G.rows();
    Eigen::Index m = CqT.cols();
    DL.setZero(num_constraints);
    if (m == 0) {
        Eigen::FullPivLU<ChMatrixDynamic<>> lu(G);
        if (!lu.isInvertible())
            return false;
        D = lu.solve(R);
        return true;
    }
    ChMatrixDynamic<> K(n + m, n + m);
    K.setZero();
    K.topLeftCorner(n, n) = G;
    K.topRightCorner(n, m) = CqT;
    K.bottomLeftCorner(m, n) = CqT.transpose();
    ChVectorDynamic<> rhs(n + m);
    rhs.head(n) = R;
    rhs.tail(m) = -Qc.head(m);
    Eigen::FullPivLU<ChMatrixDynamic<>> lu(K);
    if (!lu.isInvertible())
        return false;
    ChVectorDynamic<> sol = lu.solve(rhs);
    D = sol.head(n);
    DL.head(m) = -sol.tail(m);
    return true;
}

}  // namespace

void ChIntegrable::LoadResidual_Hv(ChVectorDynamic<>& R, const ChVectorDynamic<>& v, const double c) {
    R += c * v;
}

void ChIntegrable::LoadResidual_F(ChVectorDynamic<>& R, const double c) {
    ChState y(GetNumCoordsPosLevel() + GetNumCoordsVelLevel(), this);
    y.fill(0.0);
    double T = 0;
    StateGather(y, T);
    ChStateDelta dydt(GetNumCoordsVelLevel() + GetNumCoordsAccLevel(), this);
    ChVectorDynamic<> L(GetNumConstraints());
    if (!StateSolve(dydt, L, y, T, 0, false, UpdateFlags::UPDATE_ALL_NO_VISUAL))
        throw std::runtime_error("ChIntegrable::LoadResidual_F: StateSolve failed.");
    ChVectorDynamic<> Hdydt(dydt.size());
    Hdydt.setZero();
    LoadResidual_Hv(Hdydt, dydt, 1.0);
    R += c * Hdydt;
}

void ChIntegrable::LoadResidual_CqL(ChVectorDynamic<>& R, const ChVectorDynamic<>& L, const double c) {}

void ChIntegrable::LoadConstraint_C(ChVectorDynamic<>& Qc,
                                    const double c,
                                    const double c_vel,
                                    const bool do_clamp,
                                    const double mclam) {}

void ChIntegrable::LoadConstraint_Ct(ChVectorDynamic<>& Qc, const double c, const double c_vel) {}

bool ChIntegrable::StateSolveCorrection(ChStateDelta& Dy,
                                        ChVectorDynamic<>& DL,
                                        const ChVectorDynamic<>& R,
                                        const ChVectorDynamic<>& Qc,
                                        const double a,
                                        const double b,
                                        const ChState& y,
                                        const double T,
                                        const double dt,
                                        bool force_state_scatter,
                                        UpdateFlags update_flags,
                                        bool call_setup,
                                        bool call_analyze) {
    return SolveCorrectionFirstOrderFD(Dy, DL, R, Qc, a, b, y, T, dt, update_flags, true);
}

bool ChIntegrable::SolveCorrectionFirstOrderFD(ChStateDelta& Dy,
                                               ChVectorDynamic<>& DL,
                                               const ChVectorDynamic<>& R,
                                               const ChVectorDynamic<>& Qc,
                                               const double a,
                                               const double b,
                                               const ChState& y,
                                               const double T,
                                               const double dt,
                                               UpdateFlags update_flags,
                                               bool use_constraints) {
    const int n = (int)(GetNumCoordsVelLevel() + GetNumCoordsAccLevel());
    const int m = (int)GetNumConstraints();

    // Assemble H column by column
    ChMatrixDynamic<> H(n, n);
    ChVectorDynamic<> e(n);
    ChVectorDynamic<> col(n);
    for (int j = 0; j < n; j++) {
        e.setZero();
        e(j) = 1;
        col.setZero();
        LoadResidual_Hv(col, e, 1.0);
        H.col(j) = col;
    }

    ChMatrixDynamic<> G = a * H;

    // dF/dy = H * d(dy/dt)/dy, by forward differences of StateSolve
    if (b != 0) {
        ChStateDelta f0(n, this);
        ChStateDelta fp(n, this);
        ChVectorDynamic<> L0(m);
        ChVectorDynamic<> Lp(m);
        if (!StateSolve(f0, L0, y, T, dt, true, update_flags))
            return false;
        ChMatrixDynamic<> J(n, n);
        ChStateDelta dy(n, this);
        ChState yp(y.size(), this);
        double ymax = y.size() > 0 ? y.lpNorm<Eigen::Infinity>() : 0.0;
        for (int j = 0; j < n; j++) {
            double h = FDStep(ymax);
            dy.fill(0.0);
            dy(j) = h;
            StateIncrement(yp, y, dy);
            if (!StateSolve(fp, Lp, yp, T, dt, true, update_flags))
                return false;
            J.col(j) = (fp - f0) / h;
        }
        StateScatter(y, T, update_flags);  // Restore the original state
        G += b * H * J;
    }

    // Constraint Jacobian columns Cq' (all zero unless LoadResidual_CqL is overridden)
    ChMatrixDynamic<> CqT(n, 0);
    if (use_constraints && m > 0) {
        ChMatrixDynamic<> C(n, m);
        ChVectorDynamic<> l(m);
        for (int i = 0; i < m; i++) {
            l.setZero();
            l(i) = 1;
            col.setZero();
            LoadResidual_CqL(col, l, 1.0);
            C.col(i) = col;
        }
        if (C.cwiseAbs().maxCoeff() > 0)
            CqT = C;
    }

    ChVectorDynamic<> D;
    if (!SolveDenseKKT(G, CqT, R, Qc, D, DL, m))
        return false;
    Dy.resize(n);
    Dy = D;
    return true;
}

// -----------------------------------------------------------------------------
// Default implementations of the ChIntegrableIIorder implicit integration interface
// -----------------------------------------------------------------------------

void ChIntegrableIIorder::LoadResidual_Mv(ChVectorDynamic<>& R, const ChVectorDynamic<>& w, const double c) {
    R += c * w;
}

void ChIntegrableIIorder::LoadResidual_dFdv(ChVectorDynamic<>& R, const ChVectorDynamic<>& w, const double c) {
    const int n = (int)GetNumCoordsVelLevel();
    double wmax = w.size() > 0 ? w.lpNorm<Eigen::Infinity>() : 0.0;
    if (n == 0 || wmax == 0)
        return;

    ChState x(GetNumCoordsPosLevel(), this);
    ChStateDelta v(n, this);
    double T = 0;
    StateGather(x, v, T);

    ChVectorDynamic<> F0(n);
    F0.setZero();
    LoadResidual_F(F0, 1.0);

    // dF/dv*w = (F(x, v + h*w) - F(x, v)) / h
    double h = FDStep(v.lpNorm<Eigen::Infinity>()) / wmax;
    ChVectorDynamic<> hw = h * w;
    ChStateDelta vp = v + hw;
    StateScatter(x, vp, T, UpdateFlags::UPDATE_ALL_NO_VISUAL);
    ChVectorDynamic<> Fp(n);
    Fp.setZero();
    LoadResidual_F(Fp, 1.0);
    StateScatter(x, v, T, UpdateFlags::UPDATE_ALL_NO_VISUAL);  // Restore the original state

    R += (c / h) * (Fp - F0);
}

void ChIntegrableIIorder::LoadLumpedMass_Md(ChVectorDynamic<>& Md, double& err, const double c) {
    const int n = (int)GetNumCoordsVelLevel();
    ChVectorDynamic<> e(n);
    ChVectorDynamic<> col(n);
    double offdiag2 = 0;
    for (int j = 0; j < n; j++) {
        e.setZero();
        e(j) = 1;
        col.setZero();
        LoadResidual_Mv(col, e, 1.0);
        // Row-sum lumping (M is symmetric: column sums equal row sums)
        Md += c * col;
        offdiag2 += col.squaredNorm() - col(j) * col(j);
    }
    if (offdiag2 > 0)
        err += std::sqrt(offdiag2);
}

void ChIntegrableIIorder::LoadResidual_F(ChVectorDynamic<>& R, const double c) {
    if (m_default_F_active)
        throw std::logic_error(
            "ChIntegrableIIorder: a derived class must override either StateSolveA or LoadResidual_F.");
    m_default_F_active = true;

    const int n = (int)GetNumCoordsVelLevel();
    ChState x(GetNumCoordsPosLevel(), this);
    ChStateDelta v(n, this);
    ChStateDelta acc(n, this);
    x.fill(0.0);
    v.fill(0.0);
    double T = 0;
    StateGather(x, v, T);
    ChVectorDynamic<> L(GetNumConstraints());
    bool ok = false;
    try {
        ok = StateSolveA(acc, L, x, v, T, 0, false, UpdateFlags::UPDATE_ALL_NO_VISUAL);
    } catch (...) {
        m_default_F_active = false;
        throw;
    }
    m_default_F_active = false;
    if (!ok)
        throw std::runtime_error("ChIntegrableIIorder::LoadResidual_F: StateSolveA failed.");

    ChVectorDynamic<> Ma(n);
    Ma.setZero();
    LoadResidual_Mv(Ma, acc, 1.0);
    R += c * Ma;
}

void ChIntegrableIIorder::LoadResidual_CqL(ChVectorDynamic<>& R, const ChVectorDynamic<>& L, const double c) {}

void ChIntegrableIIorder::LoadConstraint_C(ChVectorDynamic<>& Qc,
                                           const double c,
                                           const double c_vel,
                                           const bool do_clamp,
                                           const double mclam) {}

void ChIntegrableIIorder::LoadConstraint_Ct(ChVectorDynamic<>& Qc, const double c, const double c_vel) {}

bool ChIntegrableIIorder::StateSolveCorrection(ChStateDelta& Dv,
                                               ChVectorDynamic<>& L,
                                               const ChVectorDynamic<>& R,
                                               const ChVectorDynamic<>& Qc,
                                               const double c_a,
                                               const double c_v,
                                               const double c_x,
                                               const ChState& x,
                                               const ChStateDelta& v,
                                               const double T,
                                               bool force_state_scatter,
                                               UpdateFlags update_flags,
                                               bool call_setup,
                                               bool call_analyze) {
    const int n = (int)GetNumCoordsVelLevel();
    const int m = (int)GetNumConstraints();

    if (force_state_scatter)
        StateScatter(x, v, T, update_flags);

    // Assemble M column by column
    ChMatrixDynamic<> M(n, n);
    ChVectorDynamic<> e(n);
    ChVectorDynamic<> col(n);
    for (int j = 0; j < n; j++) {
        e.setZero();
        e(j) = 1;
        col.setZero();
        LoadResidual_Mv(col, e, 1.0);
        M.col(j) = col;
    }
    ChMatrixDynamic<> H = c_a * M;

    // dF/dx and dF/dv by forward differences of LoadResidual_F
    if (c_x != 0 || c_v != 0) {
        ChVectorDynamic<> F0(n);
        ChVectorDynamic<> Fp(n);
        F0.setZero();
        LoadResidual_F(F0, 1.0);
        ChState xp(x.size(), this);
        ChStateDelta vp(v);
        ChStateDelta dx(n, this);
        double xmax = x.size() > 0 ? x.lpNorm<Eigen::Infinity>() : 0.0;
        for (int j = 0; j < n; j++) {
            if (c_x != 0) {
                double h = FDStep(xmax);
                dx.fill(0.0);
                dx(j) = h;
                StateIncrementX(xp, x, dx);
                StateScatter(xp, v, T, update_flags);
                Fp.setZero();
                LoadResidual_F(Fp, 1.0);
                H.col(j) += c_x * (Fp - F0) / h;
            }
            if (c_v != 0) {
                double h = FDStep(v(j));
                vp = v;
                vp(j) += h;
                StateScatter(x, vp, T, update_flags);
                Fp.setZero();
                LoadResidual_F(Fp, 1.0);
                H.col(j) += c_v * (Fp - F0) / h;
            }
        }
        StateScatter(x, v, T, update_flags);  // Restore the original state
    }

    // Constraint Jacobian columns Cq' (all zero unless LoadResidual_CqL is overridden)
    ChMatrixDynamic<> CqT(n, 0);
    if (m > 0) {
        LoadConstraintJacobians();
        ChMatrixDynamic<> C(n, m);
        ChVectorDynamic<> l(m);
        for (int i = 0; i < m; i++) {
            l.setZero();
            l(i) = 1;
            col.setZero();
            LoadResidual_CqL(col, l, 1.0);
            C.col(i) = col;
        }
        if (C.cwiseAbs().maxCoeff() > 0)
            CqT = C;
    }

    ChVectorDynamic<> D;
    if (!SolveDenseKKT(H, CqT, R, Qc, D, L, m))
        return false;
    Dv.resize(n);
    Dv = D;
    return true;
}

bool ChIntegrableIIorder::StateSolveCorrection(ChStateDelta& Dy,
                                               ChVectorDynamic<>& DL,
                                               const ChVectorDynamic<>& R,
                                               const ChVectorDynamic<>& Qc,
                                               const double a,
                                               const double b,
                                               const ChState& y,
                                               const double T,
                                               const double dt,
                                               bool force_state_scatter,
                                               UpdateFlags update_flags,
                                               bool call_setup,
                                               bool call_analyze) {
    // 1st order form y = {x,v}, dy/dt = {v,a}: H = identity, constraints are resolved within StateSolveA
    return SolveCorrectionFirstOrderFD(Dy, DL, R, Qc, a, b, y, T, dt, update_flags, false);
}

}  // end namespace chrono
