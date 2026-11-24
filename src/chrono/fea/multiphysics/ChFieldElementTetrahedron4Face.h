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

#ifndef CHFIELDELEMENTTETRAHEDRON4FACE_H
#define CHFIELDELEMENTTETRAHEDRON4FACE_H

#include <algorithm>

#include "chrono/core/ChQuadrature.h"
#include "chrono/fea/multiphysics/ChFieldElementTetrahedron4.h"
#include "chrono/fea/multiphysics/ChField.h"
#include "chrono/fea/multiphysics/ChFieldData.h"
#include "chrono/physics/ChLoadable.h"
#include "chrono/fea/multiphysics/ChNodeFEAfieldXYZ.h"

namespace chrono {
namespace fea {

/// @addtogroup chrono_fea
/// @{



/// Tetrahedron face. Useful for visualization etc.
/// When wrapped in a ChFieldElementLoadableSurface, this can be also used to apply UV loads 
/// of ChLoaderUVdistributed type, for example ChLoaderPressure
/// if referencing the ChFEModelDeformation, or Neumann boundary in general like imposed heat flow, etc.
/// Corner nodes, obtainable with GetNode(), are in counterclockwise order seen from the outside.

class ChApi ChFieldTetrahedron4Face : public ChFieldElementSurface {
public:
    /// Construct the specified face (0 <= id <= 3) on the given hexahedral element.
    /// The face id is: 
    /// id = 0 for face opposed to vertex 0
    /// id = 1 for face opposed to vertex 1 
    /// id = 2 for face opposed to vertex 2 
    /// id = 3 for face opposed to vertex 3
    ChFieldTetrahedron4Face(std::shared_ptr<ChFieldElementTetrahedron4> element, char id) : m_face_id(id), m_element(element) {}

    ~ChFieldTetrahedron4Face() {}

    //
    // INTERFACE to  ChFieldElement 
    //

    virtual unsigned int GetNumNodes() override { return 3; }

    virtual std::shared_ptr<ChNodeFEAbase> GetNode(unsigned int i) override {

        return m_element->GetNode(GetFaceNodeMapping(this->m_face_id)[i]);
    }

    virtual int GetSpatialDimensions() const override { return 3; }

    /// Compute shape function for the face. Note, only U,V of eta are used, assuming the two 
    /// parametric coords of the face.
    virtual void ComputeN(const ChVector3d eta, ChRowVectorDynamic<>& N) override {
        // Btw one could compute the shape functions directly as 
        //  N(0) = 1.0 - eta.x() - eta.y();
        //  N(1) = eta.x();
        //  N(2) = eta.y();
        // But here is more general, using the volume shape functions:

        ChVector3d eta_vol = ConvertFaceToVolumeCoords(m_face_id, eta);

        ChRowVectorDynamic<double> N_vol(4);
        this->m_element->ComputeN(eta_vol, N_vol);

        auto node_map = GetFaceNodeMapping(m_face_id);
        for (int i = 0; i < 3; ++i) {
            N[i] = N_vol[node_map[i]];
        }
    }
    /// Compute the parametric derivatives dN/d(u,v) of the face shape functions (2x3 matrix: rows u, v; columns face
    /// nodes). They are obtained from the volume derivatives through the chain rule of the face parametrization.
    virtual void ComputedNde(const ChVector3d eta, ChMatrixDynamic<>& dNde) override {
        ChVector3d eta_vol = ConvertFaceToVolumeCoords(m_face_id, eta);
        ChMatrixDynamic<> dNde_vol(3, 4);
        m_element->ComputedNde(eta_vol, dNde_vol);
        ChVector3d du;
        ChVector3d dv;
        GetFaceTangents(m_face_id, du, dv);
        auto node_map = GetFaceNodeMapping(m_face_id);
        dNde.setZero(2, 3);
        for (int i = 0; i < 3; ++i) {
            ChVector3d g(dNde_vol(0, node_map[i]), dNde_vol(1, node_map[i]), dNde_vol(2, node_map[i]));
            dNde(0, i) = Vdot(du, g);
            dNde(1, i) = Vdot(dv, g);
        }
    }

    /// Compute the face Jacobian. Its columns are dX/du, dX/dv and, as auxiliary third column, the unit normal
    /// (dX/du x dX/dv, normalized). Returns the area scale of the face, |dX/du x dX/dv|.
    virtual double ComputeJ(const ChVector3d eta, ChMatrix33d& J) override {
        ChVector3d eta_vol = ConvertFaceToVolumeCoords(m_face_id, eta);
        ChMatrix33d J_vol;
        m_element->ComputeJ(eta_vol, J_vol);
        ChVector3d du;
        ChVector3d dv;
        GetFaceTangents(m_face_id, du, dv);
        ChVector3d tu(J_vol * du.eigen());
        ChVector3d tv(J_vol * dv.eigen());
        ChVector3d n = Vcross(tu, tv);
        double area = n.Length();
        J.col(0) = tu.eigen();
        J.col(1) = tv.eigen();
        J.col(2) = (area > 0 ? n / area : VNULL).eigen();
        return area;
    }

    /// Compute the inverse of the face Jacobian from ComputeJ (including the normal column).
    /// Returns the determinant of the Jacobian.
    virtual double ComputeJinv(const ChVector3d eta, ChMatrix33d& Jinv) override {
        ChMatrix33d J;
        double det = ComputeJ(eta, J);
        Jinv = (det > 0) ? ChMatrix33d(J.inverse()) : ChMatrix33d(ChMatrix33d::Zero());
        return det;
    }

    /// Compute the spatial derivatives dN/dX of the face shape functions (surface gradient, 3x3 matrix with the face
    /// nodes as columns), as J^{-T} applied to the parametric derivatives. The derivative along the normal is zero.
    virtual void ComputedNdX(const ChVector3d eta, ChMatrixDynamic<>& dNdX) override {
        ChMatrixDynamic<> dNde;
        ComputedNde(eta, dNde);
        ChMatrix33d Jinv;
        ComputeJinv(eta, Jinv);
        ChMatrixDynamic<> dNde3;
        dNde3.setZero(3, 3);
        dNde3.topRows(2) = dNde;
        dNdX = Jinv.transpose() * dNde3;
    }

    /// Get the number of Gauss points of the triangle quadrature of the given order.
    virtual int GetNumQuadraturePointsForOrder(const int order) const override {
        if (order <= 1)
            return 1;
        ChQuadratureTablesTriangle* mtables = ChQuadrature::GetStaticTablesTriangle();
        int o = std::min(order, (int)mtables->Weight.size());
        return (int)mtables->Weight[o - 1].size();
    }

    /// Weight and face coordinates (u, v, 0) of the i-th Gauss point of the triangle quadrature of the given order.
    /// The weights sum to 1/2 (the parametric area), so that sum(weight * ComputeJ) is the face area.
    virtual void GetMaterialPointWeight(const int order,
                                        const int i,
                                        double& weight,
                                        ChVector3d& coords) const override {
        if (order <= 1) {
            coords.Set(CH_1_3, CH_1_3, 0.0);
            weight = 0.5;
            return;
        }
        ChQuadratureTablesTriangle* mtables = ChQuadrature::GetStaticTablesTriangle();
        int o = std::min(order, (int)mtables->Weight.size());
        coords.Set(mtables->LrootsU[o - 1][i], mtables->LrootsV[o - 1][i], 0.0);
        weight = mtables->Weight[o - 1][i] * 0.5;  // The 1/2 factor is not in the table
    }


    virtual bool IsTriangleIntegrationCompatible() const override { return true; }

    virtual ChVector3d ComputeNormal(const double U, const double V) override {
        ChVector3d p0 = *(std::dynamic_pointer_cast<ChNodeFEAfieldXYZ>(this->GetNode(0)));
        ChVector3d p1 = *(std::dynamic_pointer_cast<ChNodeFEAfieldXYZ>(this->GetNode(1)));
        ChVector3d p2 = *(std::dynamic_pointer_cast<ChNodeFEAfieldXYZ>(this->GetNode(2)));
        return Vcross(p1 - p0, p2 - p0).GetNormalized();
    }


private:

    // utility
    static std::array<int, 3>& GetFaceNodeMapping(int faceId) {
        static std::array<int, 3> ifa0{ 1, 3, 2 }; // Face opposite node 0 uses nodes 1,2,3
        static std::array<int, 3> ifa1{ 0, 2, 3 }; // Face opposite node 1 uses nodes 0,2,3
        static std::array<int, 3> ifa2{ 0, 3, 1 }; // Face opposite node 2 uses nodes 0,1,3
        static std::array<int, 3> ifa3{ 0, 1, 2 }; // Face opposite node 3 uses nodes 0,1,2
        switch (faceId) {
        case 0: return ifa0;
        case 1: return ifa1;
        case 2: return ifa2;
        case 3: return ifa3;
        default: return ifa3;
        }
    }

    // Utility: derivatives of the volume coordinates with respect to the face coordinates u, v
    static void GetFaceTangents(int faceId, ChVector3d& du, ChVector3d& dv) {
        switch (faceId) {
            case 0:  // (0, u, v)
                du.Set(0, 1, 0);
                dv.Set(0, 0, 1);
                break;
            case 1:  // (u, 0, v)
                du.Set(1, 0, 0);
                dv.Set(0, 0, 1);
                break;
            case 2:  // (u, v, 0)
                du.Set(1, 0, 0);
                dv.Set(0, 1, 0);
                break;
            default:  // (u, v, 1-u-v)
                du.Set(1, 0, -1);
                dv.Set(0, 1, -1);
                break;
        }
    }

    // utility: from  eta_face(u,v,0)  to  eta_volume(..,..,..)
    static ChVector3d ConvertFaceToVolumeCoords(int faceId, const ChVector3d& eta_face) {
        switch (faceId) {
        case 0: // Face 0: parametric coords are (eta1, eta2)
            return ChVector3d(0.0, eta_face[0], eta_face[1]);
        case 1: // Face 1: parametric coords are (eta0, eta2)
            return ChVector3d(eta_face[0], 0.0, eta_face[1]);
        case 2: // Face 2: parametric coords are (eta0, eta1)  
            return ChVector3d(eta_face[0], eta_face[1], 0.0);
        case 3: // Face 3: parametric coords are (eta0, eta1)
            return ChVector3d(eta_face[0], eta_face[1], 1.0 - eta_face[0] - eta_face[1]);
        default:
            return VNULL;
        }
    }

    char m_face_id;                                  ///< id of the face on the tetrahedron
    std::shared_ptr<ChFieldElementTetrahedron4> m_element;  ///< associated tetrahedron element
    //std::shared_ptr<ChFieldBase> m_field;  ///< associated field
};






/// @} chrono_fea

}  // end namespace fea

}  // end namespace chrono

#endif
