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
// Unit test for Bullet box-box contacts with a collision envelope: the contact
// points must lie on the faces of the boxes (not on the faces moved outward by
// the envelope) and the contact distance must be the true one, for a cube
// resting on a face, on an edge, and for two cubes touching with crossed edges.
//
// =============================================================================

#include <algorithm>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChContactContainer.h"
#include "chrono/collision/ChCollisionShapeCone.h"
#include "chrono/collision/ChCollisionShapeTriangleMesh.h"
#include "chrono/geometry/ChTriangleMeshConnected.h"

#include "gtest/gtest.h"

using namespace chrono;

class BoxContactReporter : public ChContactContainer::ReportContactCallback {
  public:
    BoxContactReporter(std::shared_ptr<ChBody> cube, double half) : m_cube(cube), m_half(half) {}

    virtual bool OnReportContact(const ChVector3d& pA,
                                 const ChVector3d& pB,
                                 const ChMatrix33<>& plane_coord,
                                 double distance,
                                 double eff_radius,
                                 const ChVector3d& react_forces,
                                 const ChVector3d& react_torques,
                                 ChContactable* contactobjA,
                                 ChContactable* contactobjB,
                                 int constraint_offset) override {
        if (contactobjA != m_cube.get() && contactobjB != m_cube.get())
            return true;
        num_contacts++;
        min_dist = std::min(min_dist, distance);
        max_dist = std::max(max_dist, distance);
        // contact point on the cube, in the cube frame: it must lie within the cube
        const ChVector3d& p = (contactobjA == m_cube.get()) ? pA : pB;
        ChVector3d loc = m_cube->TransformPointParentToLocal(p);
        double extent = std::max({std::abs(loc.x()), std::abs(loc.y()), std::abs(loc.z())});
        max_outside = std::max(max_outside, extent - m_half);
        max_lateral = std::max(max_lateral, std::max(std::abs(loc.x()), std::abs(loc.y())));
        return true;
    }

    int num_contacts = 0;
    double min_dist = 1e30;
    double max_dist = -1e30;
    double max_outside = -1e30;
    double max_lateral = 0;

  private:
    std::shared_ptr<ChBody> m_cube;
    double m_half;
};

class MeshContactReporter : public ChContactContainer::ReportContactCallback {
  public:
    MeshContactReporter(ChContactable* body) : m_body(body) {}

    virtual bool OnReportContact(const ChVector3d& pA,
                                 const ChVector3d& pB,
                                 const ChMatrix33<>& plane_coord,
                                 double distance,
                                 double eff_radius,
                                 const ChVector3d& react_forces,
                                 const ChVector3d& react_torques,
                                 ChContactable* contactobjA,
                                 ChContactable* contactobjB,
                                 int constraint_offset) override {
        if (contactobjA == m_body || contactobjB == m_body)
            num_contacts++;
        return true;
    }

    int num_contacts = 0;

  private:
    ChContactable* m_body;
};

// Cube of side L with rotation q on a fixed body A; the lowest point of the cube is 'depth' below the touching height.
static void CheckCube(std::function<std::shared_ptr<ChBody>(std::shared_ptr<ChContactMaterial>)> makeA,
                      double L,
                      const ChQuaterniond& q,
                      double z_touch,
                      double envelope,
                      const std::string& name) {
    const double depth = 1e-5;

    // the system constructor sets the default envelope and margin: set them afterwards, before creating the bodies
    ChSystemNSC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    ChCollisionModel::SetDefaultSuggestedEnvelope(envelope);
    ChCollisionModel::SetDefaultSuggestedMargin(std::min(0.01, 0.1 * L));
    auto mat = chrono_types::make_shared<ChContactMaterialNSC>();

    auto A = makeA(mat);
    A->SetFixed(true);
    sys.Add(A);

    auto cube = chrono_types::make_shared<ChBodyEasyBox>(L, L, L, 1000, false, true, mat);
    cube->SetRot(q);
    cube->SetPos(ChVector3d(0, 0, z_touch - depth));
    sys.Add(cube);

    sys.GetCollisionSystem()->BindAll();
    sys.ComputeCollisions();

    auto reporter = chrono_types::make_shared<BoxContactReporter>(cube, L / 2);
    sys.GetContactContainer()->ReportAllContacts(reporter);

    std::cout << name << "  L " << L << "  envelope " << envelope << "  contacts " << reporter->num_contacts << "  min distance " << reporter->min_dist
              << "  max distance of a point outside the cube " << reporter->max_outside << std::endl;

    ASSERT_GE(reporter->num_contacts, 1);
    ASSERT_NEAR(reporter->min_dist, -depth, 0.05 * L * 1e-3 + 2 * depth);
    ASSERT_LT(reporter->max_outside, 2e-7 + 1e-6 * L);  // rounding errors of single precision Bullet
}

TEST(BulletCollision, BoxBoxEnvelope) {
    using Mat = std::shared_ptr<ChContactMaterial>;
    for (double envelope : {0.03, 0.005, 0.0}) {
        for (double L : {0.05, 0.5}) {
            // a floor box with its top face at z = 0
            auto floor = [](Mat mat) {
                auto body = chrono_types::make_shared<ChBodyEasyBox>(4.0, 4.0, 0.2, 1000, false, true, mat);
                body->SetPos(ChVector3d(0.3, -0.2, -0.1));
                return body;
            };
            // cube resting on a face, yawed
            for (double yaw : {0.0, 30.0})
                CheckCube(floor, L, QuatFromAngleZ(yaw * CH_DEG_TO_RAD), L / 2, envelope, "cube on a face, yaw " + std::to_string(yaw));
            // cube resting on an edge
            for (double roll : {10.0, 45.0}) {
                double t = roll * CH_DEG_TO_RAD;
                CheckCube(floor, L, QuatFromAngleX(t), L / 2 * (std::cos(t) + std::sin(t)), envelope, "cube on an edge, roll " + std::to_string(roll));
            }
            // crossed edges: lower cube rotated 45 deg about Y, upper cube 45 deg about X
            auto lower = [L](Mat mat) {
                auto body = chrono_types::make_shared<ChBodyEasyBox>(L, L, L, 1000, false, true, mat);
                body->SetRot(QuatFromAngleY(CH_PI_4));
                return body;
            };
            CheckCube(lower, L, QuatFromAngleX(CH_PI_4), L * std::sqrt(2.0), envelope, "crossed edges");
        }
    }
}

// A flat face resting on a flat face must be supported over the face: a convex hull cube, a cylinder on its end and a
// cone on its base, resting on a floor box, get at least three contact points at the true depth, spread over the
// face up to the edge rounding by the safe margin.
static void CheckFlatSupport(const std::string& type, double L, double envelope, double margin) {
    const double depth = 1e-5;

    ChSystemNSC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    ChCollisionModel::SetDefaultSuggestedEnvelope(envelope);
    ChCollisionModel::SetDefaultSuggestedMargin(margin);
    auto mat = chrono_types::make_shared<ChContactMaterialNSC>();

    auto floor = chrono_types::make_shared<ChBodyEasyBox>(4.0, 4.0, 0.2, 1000, false, true, mat);
    floor->SetPos(ChVector3d(0.3, -0.2, -0.1));
    floor->SetFixed(true);
    sys.Add(floor);

    std::shared_ptr<ChBody> body;
    double z = L / 2;  // center height when touching
    if (type == "hull") {
        std::vector<ChVector3d> points;
        for (int i = 0; i < 8; i++)
            points.push_back(ChVector3d(i & 1 ? L / 2 : -L / 2, i & 2 ? L / 2 : -L / 2, i & 4 ? L / 2 : -L / 2));
        body = chrono_types::make_shared<ChBodyEasyConvexHull>(points, 1000, false, true, mat);
    } else if (type == "cylinder") {
        body = chrono_types::make_shared<ChBodyEasyCylinder>(ChAxis::Z, L / 2, L, 1000, false, true, mat);
    } else {
        body = chrono_types::make_shared<ChBody>();
        body->AddCollisionShape(chrono_types::make_shared<ChCollisionShapeCone>(mat, L / 2, L));
        body->EnableCollision(true);
    }
    body->SetPos(ChVector3d(0, 0, z - depth));
    sys.Add(body);

    sys.GetCollisionSystem()->BindAll();
    for (int i = 0; i < 3; i++)
        sys.ComputeCollisions();

    auto reporter = chrono_types::make_shared<BoxContactReporter>(body, L / 2);
    sys.GetContactContainer()->ReportAllContacts(reporter);

    std::cout << type << "  L " << L << "  envelope " << envelope << "  margin " << margin << "  contacts " << reporter->num_contacts << "  distances [" << reporter->min_dist
              << ", " << reporter->max_dist << "]  lateral extent " << reporter->max_lateral << std::endl;

    ASSERT_GE(reporter->num_contacts, 3);
    ASSERT_NEAR(reporter->min_dist, -depth, 1e-6 + 1e-6 * L);
    ASSERT_NEAR(reporter->max_dist, -depth, 1e-6 + 1e-6 * L);
    // points within the body and out to its edges, less the rounding of the edges by the safe margin
    ASSERT_LT(reporter->max_outside, 2e-7 + 1e-6 * L);  // rounding errors of single precision Bullet
    ASSERT_GT(reporter->max_lateral, L / 2 - 2 * margin);
}

TEST(BulletCollision, FlatFaceSupport) {
    for (std::string type : {"hull", "cylinder", "cone"})
        for (double L : {0.05, 0.5})
            for (double envelope : {0.03, 0.005, 0.0})
                CheckFlatSupport(type, L, envelope, std::min(0.01, 0.02 * L));
}

// Flat square mesh of side 'size' in the plane z = 0, with n x n cells of two triangles each
static std::shared_ptr<ChTriangleMeshConnected> FlatMesh(double size, int n) {
    auto mesh = chrono_types::make_shared<ChTriangleMeshConnected>();
    for (int j = 0; j <= n; j++)
        for (int i = 0; i <= n; i++)
            mesh->m_vertices.push_back(ChVector3d(-size / 2 + i * size / n, -size / 2 + j * size / n, 0));
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            int a = j * (n + 1) + i;
            mesh->m_face_v_indices.push_back(ChVector3i(a, a + 1, a + n + 2));
            mesh->m_face_v_indices.push_back(ChVector3i(a, a + n + 2, a + n + 1));
        }
    return mesh;
}

// The triangles of a mesh are elements of a discretised surface, supported at their corners by the neighbouring
// triangles: a flat mesh of small triangles resting on a floor box gets one contact point per triangle. A box resting
// on a single large mesh triangle is still supported over its face (by at least three points).
TEST(BulletCollision, MeshTriangleSupport) {
    const double depth = 1e-5;
    const int n = 4;

    for (double envelope : {0.03, 0.005}) {
        ChSystemNSC sys;
        sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
        ChCollisionModel::SetDefaultSuggestedEnvelope(envelope);
        ChCollisionModel::SetDefaultSuggestedMargin(0.002);
        auto mat = chrono_types::make_shared<ChContactMaterialNSC>();

        // small triangles (mesh 0.2 x 0.2) resting on a floor box with its top face at z = 0
        auto floor = chrono_types::make_shared<ChBodyEasyBox>(4.0, 4.0, 0.2, 1000, false, true, mat);
        floor->SetPos(ChVector3d(0.3, -0.2, -0.1));
        floor->SetFixed(true);
        sys.Add(floor);
        auto plate = chrono_types::make_shared<ChBody>();
        plate->AddCollisionShape(chrono_types::make_shared<ChCollisionShapeTriangleMesh>(mat, FlatMesh(0.2, n), false, false, 0.0));
        plate->EnableCollision(true);
        plate->SetPos(ChVector3d(0, 0, -depth));
        sys.Add(plate);

        // box of side 0.2 resting on a mesh floor of two triangles (4 x 4), inside one of them
        auto mesh_floor = chrono_types::make_shared<ChBody>();
        mesh_floor->AddCollisionShape(chrono_types::make_shared<ChCollisionShapeTriangleMesh>(mat, FlatMesh(4.0, 1), true, false, 0.0));
        mesh_floor->EnableCollision(true);
        mesh_floor->SetPos(ChVector3d(0, 0, -10));
        mesh_floor->SetFixed(true);
        sys.Add(mesh_floor);
        auto cube = chrono_types::make_shared<ChBodyEasyBox>(0.2, 0.2, 0.2, 1000, false, true, mat);
        cube->SetPos(ChVector3d(0.8, -0.8, -10 + 0.1 - depth));
        sys.Add(cube);

        sys.GetCollisionSystem()->BindAll();
        for (int i = 0; i < 3; i++)
            sys.ComputeCollisions();

        auto plate_reporter = chrono_types::make_shared<MeshContactReporter>(plate.get());
        auto cube_reporter = chrono_types::make_shared<BoxContactReporter>(cube, 0.1);
        sys.GetContactContainer()->ReportAllContacts(plate_reporter);
        sys.GetContactContainer()->ReportAllContacts(cube_reporter);

        std::cout << "envelope " << envelope << "  mesh triangles " << 2 * n * n << "  contacts " << plate_reporter->num_contacts << "  box on a mesh triangle: contacts "
                  << cube_reporter->num_contacts << "  distances [" << cube_reporter->min_dist << ", " << cube_reporter->max_dist << "]" << std::endl;

        ASSERT_EQ(plate_reporter->num_contacts, 2 * n * n);
        ASSERT_GE(cube_reporter->num_contacts, 3);
        ASSERT_NEAR(cube_reporter->min_dist, -depth, 1e-6);
        ASSERT_NEAR(cube_reporter->max_dist, -depth, 1e-6);
        ASSERT_LT(cube_reporter->max_outside, 1e-6);
    }
}
