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

#include <cmath>
#include <limits>

#include "chrono/solver/ChSolverBB.h"
#include "chrono/utils/ChConstants.h"

namespace chrono {

// Register into the object factory, to enable run-time dynamic creation and persistence
CH_FACTORY_REGISTER(ChSolverBB)

ChSolverBB::ChSolverBB() : n_armijo(10), max_armijo_backtrace(3), lastgoodres(0) {}

double ChSolverBB::Solve(ChSystemDescriptor& sysd) {
    if (!sysd.SupportsSchurComplement()) {
        std::cerr << "\n\nChSolverBB: Can NOT use Barzilai-Borwein solver if\n"
                  << " - there are stiffness or damping matrices, or\n "
                  << " - no inverse mass matrix was provided" << std::endl;
        throw std::runtime_error("ChSolverBB: System descriptor does not support Schur complement-based solvers.");
    }

    // Tuning of the spectral gradient search
    double a_min = 1e-13;
    double a_max = 1e13;
    double sigma_min = 0.1;
    double sigma_max = 0.9;
    double alpha = 0.0001;
    double gamma = 1e-4;
    double gdiff = 0.000001;

    bool do_BB1e2 = true;
    bool do_BB1 = false;
    bool do_BB2 = false;
    double neg_BB1_fallback = 0.11;
    double neg_BB2_fallback = 0.12;

    m_iterations = 0;

    int nc = sysd.CountActiveConstraints();
    if (verbose)
        std::cout << "\n-----Barzilai-Borwein, solving nc=" << nc << "unknowns" << std::endl;

    // Allocate auxiliary vectors
    ChVectorDynamic<> ml(nc);
    ChVectorDynamic<> ml_candidate(nc);
    ChVectorDynamic<> mg(nc);
    ChVectorDynamic<> mg_p(nc);
    ChVectorDynamic<> ml_p(nc);
    ChVectorDynamic<> mdir(nc);
    ChVectorDynamic<> mb(nc);
    ChVectorDynamic<> mb_tmp(nc);
    ChVectorDynamic<> ms(nc);
    ChVectorDynamic<> my(nc);
    ChVectorDynamic<> mDg(nc);
    ChVectorDynamic<> ml_proj(nc);

    // Update auxiliary data in constraints
    // Average entries for friction constraints
    sysd.SchurComplementUpdateConstraints(true);

    // Calculate the Schur complement transformed RHS
    // Cache M^-1 * f in 'Mif'
    ChVectorDynamic<> Mif;
    sysd.SchurComplementRHS(mb, &Mif);

    // Initialize lambdas
    if (m_warm_start)
        sysd.FromConstraintsToVector(ml);
    else
        ml.setZero();

    // Initial projection of ml   ***TO DO***?
    sysd.ConstraintsProject(ml);

    // Cache the inverse diagonal of N (used by the diagonal preconditioner).
    // A row with a zero, negligible or non-finite diagonal entry acts on no active variable (e.g. a joint between two
    // fixed bodies): its multiplier is kept at zero and its right-hand side is ignored, with a unit preconditioner entry.
    ChVectorDynamic<> iD(nc);
    sysd.SchurComplementProduct(mg, ml, &iD);
    double Dmax = 0;
    for (int i = 0; i < nc; i++) {
        if (std::isfinite(1 / iD(i)))
            Dmax = std::max(Dmax, 1 / iD(i));
    }
    for (int i = 0; i < nc; i++) {
        double Dii = 1 / iD(i);
        if (!std::isfinite(Dii) || !(Dii > std::numeric_limits<double>::epsilon() * Dmax)) {
            ml(i) = 0;
            mb(i) = 0;
            iD(i) = 1;
        }
    }

    // Fallback solution
    lastgoodres = 1e30;
    ml_candidate = ml;

    // Calculate g = grad(0.5*l'*N*l-l'*b) = N*l-b

    // 1) g = N*l
    sysd.SchurComplementProduct(mg, ml);

    // 2) g = N*l -b
    mg -= mb;

    mg_p = mg;

    // initial norm of the gradient
    double mg_p_init_norm = std::max(1e-10, mg_p.norm());

    // Iterations

    // Objective at the initial iterate, f = 0.5*l'*N*l - l'*b = 0.5*l'*(g - b)
    double mf = 0.5 * ml.dot(mg - mb);
    double mf_p = mf;

    // Objective at the accepted iterates, for the non-monotone line search
    std::vector<double> f_hist;
    f_hist.push_back(mf);

    std::fill(violation_history.begin(), violation_history.end(), 0.0);
    std::fill(dlambda_history.begin(), dlambda_history.end(), 0.0);

    for (int iter = 0; iter < m_max_iterations; iter++) {
        // Dg = Di*g;
        mDg = mg;
        if (m_use_precond)
            mDg = mDg.array() * iD.array();

        // dir  = [P(l - alpha*Dg) - l]
        mdir = ml - alpha * mDg;        // dir = l - alpha*Dg
        sysd.ConstraintsProject(mdir);  // dir = P(l - alpha*Dg)
        mdir -= ml;                     // dir = P(l - alpha*Dg) - l

        // dTg = dir'*g;
        double dTg = mdir.dot(mg);

        // BB dir backward!? fallback to nonpreconditioned dir
        if (dTg > 1e-8) {
            // dir  = [P(l - alpha*g) - l]
            mdir = ml - alpha * mg;         // dir = l - alpha*g
            sysd.ConstraintsProject(mdir);  // dir = P(l - alpha*g) ...
            mdir -= ml;                     // dir = P(l - alpha*g) - l
            // dTg = d'*g;
            dTg = mdir.dot(mg);
        }

        double lambda = 1;

        int n_backtracks = 0;
        bool armijo_repeat = true;

        // Reference value of the non-monotone Armijo condition: the largest objective over the last n_armijo
        // accepted iterates, the current one included
        double max_compare = mf;
        for (int h = 1; h < std::min((int)f_hist.size(), this->n_armijo); h++)
            max_compare = std::max(max_compare, f_hist[f_hist.size() - 1 - h]);

        while (armijo_repeat) {
            // l_p = l + lambda*dir;
            ml_p = ml + lambda * mdir;

            // m_tmp = Nl_p = N*l_p;
            sysd.SchurComplementProduct(mb_tmp, ml_p);

            // g_p = N * l_p - b  = Nl_p - b
            mg_p = mb_tmp - mb;

            // f_p = 0.5*l_p'*N*l_p - l_p'*b  = l_p'*(0.5*Nl_p - b);
            mf_p = ml_p.dot(0.5 * mb_tmp - mb);

            if (mf_p > max_compare + gamma * lambda * dTg) {
                armijo_repeat = true;
                double lambdanew = -lambda * lambda * dTg / (2 * (mf_p - mf - lambda * dTg));
                lambda = std::max(sigma_min * lambda, std::min(sigma_max * lambda, lambdanew));
                if (verbose)
                    std::cout << " Repeat Armijo, new lambda=" << lambda << std::endl;
            } else {
                armijo_repeat = false;
            }

            n_backtracks = n_backtracks + 1;
            if (n_backtracks > this->max_armijo_backtrace)
                armijo_repeat = false;
        }

        ms = ml_p - ml;  // s = l_p - l;
        my = mg_p - mg;  // y = g_p - g;
        ml = ml_p;       // l = l_p;
        mg = mg_p;       // g = g_p;
        mf = mf_p;       // f = f_p;
        f_hist.push_back(mf);

        if (((do_BB1e2) && (iter % 2 == 0)) || do_BB1) {
            if (m_use_precond)
                mb_tmp = ms.array() / iD.array();
            else
                mb_tmp = ms;
            double sDs = ms.dot(mb_tmp);
            double sy = ms.dot(my);
            if (sy <= 0) {
                alpha = neg_BB1_fallback;
            } else {
                double alph = sDs / sy;  // (s,Ds)/(s,y)   BB1
                alpha = std::min(a_max, std::max(a_min, alph));
            }
        }

        /*
        // this is a modified rayleight quotient - looks like it works anyway...
        if (((do_BB1e2) && (iter%2 ==0)) || do_BB1)
        {
            double ss = ms.MatrDot(ms,ms);
            mb_tmp = my;
            if (m_use_precond)
                mb_tmp.MatrDivScale(mD);
            double sDy = ms.MatrDot(ms, mb_tmp);
            if (sDy <= 0)
            {
                alpha = neg_BB1_fallback;
            }
            else
            {
                double alph = ss / sDy;  // (s,s)/(s,Di*y)   BB1 (modified version)
                alpha = std::min (a_max, std::max(a_min, alph));
            }
        }
        */

        if (((do_BB1e2) && (iter % 2 != 0)) || do_BB2) {
            double sy = ms.dot(my);
            if (m_use_precond)
                mb_tmp = my.array() * iD.array();
            else
                mb_tmp = my;
            double yDy = my.dot(mb_tmp);
            if (sy <= 0) {
                alpha = neg_BB2_fallback;
            } else {
                double alph = sy / yDy;  // (s,y)/(y,Di*y)   BB2
                alpha = std::min(a_max, std::max(a_min, alph));
            }
        }

        // Project the gradient (for rollback strategy)
        // g_proj = (l-project_orthogonal(l - gdiff*g, fric))/gdiff;
        // evaluated as g + [(l - gdiff*g) - ProjectionOperator(l - gdiff*g)] / gdiff, exact for the entries not changed
        // by the projection (the difference l - (l - gdiff*g) loses gdiff*g in floating point when |l| >> gdiff*|g|)
        mb_tmp = ml - gdiff * mg;
        ml_proj = mb_tmp;
        sysd.ConstraintsProject(ml_proj);          // ml_proj = ProjectionOperator(l - gdiff * g)
        mb_tmp = mg + (mb_tmp - ml_proj) / gdiff;  // mb_tmp = [l - ProjectionOperator(l - gdiff * g)] / gdiff
        double g_proj_norm = mb_tmp.norm();        // infinity norm is faster..

        // Stop on a non-finite iterate, keeping the best solution found so far
        if (!std::isfinite(g_proj_norm))
            break;

        // Rollback solution: the last best candidate ('l' with lowest projected gradient)
        // in fact the method is not monotone and it is quite 'noisy', if we do not
        // do this, a prematurely truncated iteration might give a crazy result.
        if (g_proj_norm < lastgoodres) {
            lastgoodres = g_proj_norm;
            ml_candidate = ml;
        }

        // METRICS - convergence, plots, etc

        double maxdeltalambda = ms.lpNorm<Eigen::Infinity>();
        double maxd = lastgoodres / mg_p_init_norm;

        // For recording into correction/residuals/violation history, if debugging
        if (this->record_violation_history)
            AtIterationEnd(maxd, maxdeltalambda, iter);

        if (verbose)
            std::cout << "  iter=" << iter << "   f=" << mf_p << "  |d|=" << maxd << "  |s|=" << maxdeltalambda
                      << std::endl;

        m_iterations++;

        if (maxd < m_tolerance) {
            if (verbose)
                std::cout << "Converged at iter: " << iter << " with residual: " << maxd << std::endl;
            break;
        }
    }

    // Fallback to best found solution (might be useful because of nonmonotonicity)
    ml = ml_candidate;

    // Resulting DUAL variables:
    // store ml temporary vector into ChConstraint 'l_i' multipliers
    sysd.FromVectorToConstraints(ml);

    // Resulting PRIMAL variables:
    // compute the primal variables as   v = (M^-1)(f + Cq'*l)
    sysd.SchurComplementIncrementVariables(&Mif);

    if (verbose)
        std::cout << "-----" << std::endl;

    return lastgoodres;
}

// -----------------------------------------------------------------------------

void ChSolverBB::ArchiveOut(ChArchiveOut& archive_out) {
    // version number
    archive_out.VersionWrite<ChSolverBB>();
    // serialize parent class
    ChIterativeSolverVI::ArchiveOut(archive_out);
    // serialize all member data:
    archive_out << CHNVP(n_armijo);
    archive_out << CHNVP(max_armijo_backtrace);
    archive_out << CHNVP(m_use_precond);
}

void ChSolverBB::ArchiveIn(ChArchiveIn& archive_in) {
    // version number
    /*int version =*/archive_in.VersionRead<ChSolverBB>();
    // deserialize parent class
    ChIterativeSolverVI::ArchiveIn(archive_in);
    // stream in all member data:
    archive_in >> CHNVP(n_armijo);
    archive_in >> CHNVP(max_armijo_backtrace);
    archive_in >> CHNVP(m_use_precond);
}

}  // end namespace chrono
