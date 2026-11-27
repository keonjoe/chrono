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
// Authors: Alessandro Tasora 
// =============================================================================


#include "chrono/fea/multiphysics/ChLoaderHeatRadiation.h"
#include "chrono/fea/multiphysics/ChFieldElementLoadableSurface.h"

namespace chrono {
namespace fea {


void ChLoaderHeatRadiation::ComputeF(double U, double V, ChVectorDynamic<>& F, ChVectorDynamic<>* state_x, ChVectorDynamic<>* state_w) {
    F(0) = 0; // fallback if no face
    if (auto loface = std::dynamic_pointer_cast<ChFieldElementLoadableSurface>(this->loadable)) {
        if (auto face = std::dynamic_pointer_cast<ChFieldElementSurface>(loface->GetElement())) {

            // Ok, we are loading a "face" finite element, so we can use NodeData(node) and fetch associated temperature .
            // Compute shape functions for interpolating temperature from the m_temp auxiliary field:
            ChRowVectorDynamic<> N(face->GetNumNodes());
            face->ComputeN(ChVector3d(U, V, 0), N);
            // Interpolate T as  T = [T1,T2,T3,...] * N^T
            ChRowVectorDynamic<> T_hat(face->GetNumNodes());
            // Use the temperatures in state_w if provided (as when computing the Jacobian), the nodal ones otherwise.
            for (unsigned int in = 0; in < face->GetNumNodes(); ++in) {
                T_hat(in) = state_w ? (*state_w)(in) : m_temp->NodeData(face->GetNode(in)).T();
            }
            double T = T_hat * N.transpose();

            // Compute the net heat flux exchanged with the environment using the Stefan-Boltzmann law (grey surface in
            // a large enclosure, absolute temperatures): W_per_area = e * sigma * (T^4 - T_env^4), outgoing if T > T_env
            double sigma = 5.670374419e-8; // Stefan-Boltzmann constant
            F(0) = -this->m_emissivity * sigma * (std::pow(T, 4) - std::pow(m_T_env, 4));
        }
    } 
}

}  // end namespace fea
}  // end namespace chrono
