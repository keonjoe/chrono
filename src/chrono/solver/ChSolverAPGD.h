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

#ifndef CHSOLVERAPGD_H
#define CHSOLVERAPGD_H

#include "chrono/solver/ChIterativeSolverVI.h"

namespace chrono {

/// @addtogroup chrono_solver
/// @{

/// An iterative solver based on Nesterov's Projected Gradient Descent.
///
/// Notes:
/// - can not solve problems including stiffness or damping blocks in system descriptor (e.g. coming from FEA)
/// - diagonal preconditioning is available (enabled by default)
///
/// See ChSystemDescriptor for more information about the problem formulation and solver data structures.
class ChApi ChSolverAPGD : public ChIterativeSolverVI {
  public:
    ChSolverAPGD();

    ~ChSolverAPGD() {}

    virtual Type GetType() const override { return Type::APGD; }

    /// Performs the solution of the problem.
    virtual double Solve(ChSystemDescriptor& sysd) override;

    /// Return the tolerance error reached during the last solve.
    /// For the APGD solver, this is the norm of the projected gradient.
    virtual double GetError() const override { return residual; }

  private:
    double residual;
    int nc;
    ChVectorDynamic<> gamma_hat, gammaNew, g, y, gamma, yNew, r, tmp;
};

/// @} chrono_solver

}  // end namespace chrono

#endif
