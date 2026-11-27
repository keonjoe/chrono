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

#include "chrono/solver/ChSolverAPGD.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <valarray>
#include <vector>

namespace chrono {

// Register into the object factory, to enable run-time dynamic creation and persistence
CH_FACTORY_REGISTER(ChSolverAPGD)

ChSolverAPGD::ChSolverAPGD() : nc(0), residual(0.0) {}

double ChSolverAPGD::Solve(ChSystemDescriptor& sysd) {
    if (!sysd.SupportsSchurComplement()) {
        std::cerr << "\n\nChSolverAPGD: Can NOT use APGD solver if\n"
                  << " - there are stiffness or damping matrices, or\n "
                  << " - no inverse mass matrix was provided" << std::endl;
        throw std::runtime_error("ChSolverAPGD: System descriptor does not support Schur complement-based solvers.");
    }

    // Diagonal preconditioning: the iteration works in the metric of D = diag(N). The Euclidean projection onto the
    // friction cones is a projection in this metric only if D is uniform within each cone, so the diagonal entries of
    // friction constraints are averaged (not available with a non-block diagonal mass matrix).
    bool use_precond = m_use_precond && !sysd.HasMassInverse();

    // Update auxiliary data in all constraints
    sysd.SchurComplementUpdateConstraints(use_precond);

    double L, t;
    double theta;
    double thetaNew;
    double Beta;

    nc = sysd.CountActiveConstraints();
    gamma_hat.resize(nc);
    gammaNew.resize(nc);
    g.resize(nc);
    y.resize(nc);
    gamma.resize(nc);
    yNew.resize(nc);
    r.resize(nc);
    tmp.resize(nc);

    residual = 10e30;

    Beta = 0.0;

    // Compute the RHS vector in the Schur complement equation N*l = b_schur
    // Cache M^-1 * f in 'Mif'
    ChVectorDynamic<> Mif;
    sysd.SchurComplementRHS(r, &Mif);

    // If no constraints, return now. Variables contain M^-1 * f after call to SchurComplementRHS.
    // This early exit is needed, else we get division by zero and a potential infinite loop.
    if (nc == 0) {
        residual = 0;
        return 0;
    }

    // (1) gamma_0 = zeros(nc,1)
    if (m_warm_start) {
        sysd.SchurComplementIncrementVariables();
    } else {
        for (auto& constraint : sysd.GetConstraints())
            constraint->SetLagrangeMultiplier(0);
    }
    sysd.FromConstraintsToVector(gamma);

    // Inverse of the diagonal preconditioner, D^-1 = diag(N)^-1 (identity if no preconditioning).
    // A row with a zero, negligible or non-finite diagonal entry acts on no active variable (e.g. a joint between two
    // fixed bodies): its multiplier is kept at zero and its right-hand side is ignored, with a unit preconditioner entry.
    ChVectorDynamic<> Di(nc);
    sysd.SchurComplementProduct(yNew, gamma, &Di);
    double Dmax = 0;
    for (int i = 0; i < nc; i++) {
        if (std::isfinite(1 / Di(i)))
            Dmax = std::max(Dmax, 1 / Di(i));
    }
    for (int i = 0; i < nc; i++) {
        double Dii = 1 / Di(i);
        if (!std::isfinite(Dii) || !(Dii > std::numeric_limits<double>::epsilon() * Dmax)) {
            gamma(i) = 0;
            r(i) = 0;
            Di(i) = 1;
        }
    }
    if (!use_precond)
        Di.setOnes();
    ChVectorDynamic<> gstep(nc);
    ChVectorDynamic<> gproj(nc);

    // (2) gamma_hat_0 = ones(nc,1)
    gamma_hat.setConstant(1.0);

    // (3) y_0 = gamma_0
    y = gamma;

    // (4) theta_0 = 1
    theta = 1.0;

    // (5) L_k = norm(N * (gamma_0 - gamma_hat_0)) / norm(gamma_0 - gamma_hat_0)
    //     in the metric of D: norm(D^(-1/2) * N * (gamma_0 - gamma_hat_0)) / norm(D^(1/2) * (gamma_0 - gamma_hat_0))
    tmp = gamma - gamma_hat;
    sysd.SchurComplementProduct(yNew, tmp);  // yNew = N * tmp = N * (gamma - gamma_hat)
    L = std::sqrt(yNew.dot(Di.cwiseProduct(yNew)) / tmp.dot(tmp.cwiseQuotient(Di)));
    if (!std::isfinite(L) || L <= 0)
        L = 1.0;     // no estimate (e.g. N * tmp = 0); the line search adjusts L
    yNew.setZero();  //// RADU  is this really necessary here?

    // (6) t_k = 1 / L_k
    t = 1.0 / L;

    // Fallback solution, if no iterate gets a finite residual
    gamma_hat = gamma;

    //// RADU
    //// Check consistency (e.g. sign of 'r' in comments vs. code)

    std::fill(violation_history.begin(), violation_history.end(), 0.0);
    std::fill(dlambda_history.begin(), dlambda_history.end(), 0.0);

    // (7) for k := 0 to N_max
    for (m_iterations = 0; m_iterations < m_max_iterations; m_iterations++) {
        // (8) g = N * y_k - r
        // (9) gamma_(k+1) = ProjectionOperator(y_k - t_k * D^-1 * g)
        sysd.SchurComplementProduct(g, y);  // g = N * y
        g -= r;                             // g = N * y - r
        gammaNew = y - t * Di.cwiseProduct(g);
        sysd.ConstraintsProject(gammaNew);

        // (10) while 0.5 * gamma_(k+1)' * N * gamma_(k+1) - gamma_(k+1)' * r >
        //            0.5 * y_k' * N * y_k - y_k' * r + g' * (gamma_(k+1) - y_k) + 0.5 * L_k * |gamma_(k+1) - y_k|_D^2
        // For the quadratic objective, with d = gamma_(k+1) - y_k, this is d' * N * d > L_k * d' * D * d.
        // Evaluate it in this form: the difference of the objective values loses all significant digits once the
        // step is small compared to the multipliers, and the test would then double L_k without bound.
        tmp = gammaNew - y;
        sysd.SchurComplementProduct(yNew, tmp);  // yNew = N * (gamma_(k+1) - y_k)

        while (tmp.dot(yNew) > L * tmp.dot(tmp.cwiseQuotient(Di))) {
            // (11) L_k = 2 * L_k
            L = 2.0 * L;

            // (12) t_k = 1 / L_k
            t = 1.0 / L;

            // (13) gamma_(k+1) = ProjectionOperator(y_k - t_k * D^-1 * g)
            gammaNew = y - t * Di.cwiseProduct(g);
            sysd.ConstraintsProject(gammaNew);

            tmp = gammaNew - y;
            sysd.SchurComplementProduct(yNew, tmp);  // yNew = N * (gamma_(k+1) - y_k)
        }  // (14) endwhile

        // (15) theta_(k+1) = (-theta_k^2 + theta_k * sqrt(theta_k^2 + 4)) / 2
        thetaNew = (-theta * theta + theta * std::sqrt(theta * theta + 4.0)) / 2.0;

        // (16) Beta_(k+1) = theta_k * (1 - theta_k) / (theta_k^2 + theta_(k+1))
        Beta = theta * (1.0 - theta) / (theta * theta + thetaNew);

        // (17) y_(k+1) = gamma_(k+1) + Beta_(k+1) * (gamma_(k+1) - gamma_k)
        yNew = gammaNew + Beta * (gammaNew - gamma);

        // (18) r = r(gamma_(k+1))
        // Project the gradient (for rollback strategy)
        // g_proj = (l-project_orthogonal(l - gdiff*g, fric))/gdiff;
        // evaluated as g + [(l - gdiff*g) - ProjectionOperator(l - gdiff*g)] / gdiff, exact for the entries not changed
        // by the projection (the difference l - (l - gdiff*g) loses gdiff*g in floating point when |l| >> gdiff*|g|)
        double gdiff = 1.0 / (nc * nc);
        sysd.SchurComplementProduct(tmp, gammaNew);  // tmp = N * gammaNew
        tmp -= r;                                    // tmp = g = N * gammaNew - r
        gstep = gammaNew - gdiff * tmp;              // gstep = gammaNew - gdiff * g
        gproj = gstep;
        sysd.ConstraintsProject(gproj);  // gproj = ProjectionOperator(gammaNew - gdiff * g)
        double res = (tmp + (gstep - gproj) / gdiff).norm();

        // Stop on a non-finite iterate, keeping the best solution found so far
        if (!std::isfinite(res))
            break;

        if (res < residual) {      // (19) if r < epsilon_min
            residual = res;        // (20) r_min = r
            gamma_hat = gammaNew;  // (21) gamma_hat = gamma_(k+1)
        }                          // (22) endif

        if (residual < m_tolerance) {  // (23) if r < Tau
            m_iterations++;            //      (count the current iteration)
            break;                     // (24) break
        }                              // (25) endif

        if (g.dot(gammaNew - gamma) > 0) {  // (26) if g' * (gamma_(k+1) - gamma_k) > 0
            yNew = gammaNew;                // (27) y_(k+1) = gamma_(k+1)
            thetaNew = 1.0;                 // (28) theta_(k+1) = 1
        }                                   // (29) endif

        // (30) L_k = 0.9 * L_k
        L = 0.9 * L;

        // (31) t_k = 1 / L_k
        t = 1.0 / L;

        // perform some tasks at the end of the iteration
        if (this->record_violation_history) {
            AtIterationEnd(residual, (gammaNew - gamma).lpNorm<Eigen::Infinity>(), m_iterations);
        }

        // Update iterates
        theta = thetaNew;
        gamma = gammaNew;
        y = yNew;
    }  // (32) endfor

    if (verbose)
        std::cout << "Residual: " << residual << ", Iter: " << m_iterations << std::endl;

    // (33) return Value at time step t_(l+1), gamma_(l+1) := gamma_hat
    sysd.FromVectorToConstraints(gamma_hat);

    // Resulting PRIMAL variables:
    // compute the primal variables as   v = (M^-1)(f + Cq'*l)
    sysd.SchurComplementIncrementVariables(&Mif);

    return residual;
}

}  // end namespace chrono
