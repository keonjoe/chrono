// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2026 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
//
// Tests for mesh files that cannot be loaded and for mesh visualization of the Marder tracked vehicle:
// - a null visualization mesh is ignored
// - a visualization mesh file that cannot be loaded falls back to the primitive shapes
// - a collision mesh that cannot be loaded is skipped
// - a missing terrain or tire contact mesh file throws an exception that names the file
// - the Marder vehicle can be created and initialized with mesh visualization on all of its parts
//
// =============================================================================

#include <filesystem>
#include <fstream>
#include <string>

#include "gtest/gtest.h"

#include "chrono/assets/ChVisualShapeTriangleMesh.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/utils/ChBodyGeometry.h"

#include "chrono_vehicle/ChVehicleDataPath.h"
#include "chrono_vehicle/terrain/RigidTerrain.h"
#include "chrono_vehicle/terrain/SCMTerrain.h"
#include "chrono_vehicle/wheeled_vehicle/ChWheel.h"
#include "chrono_vehicle/wheeled_vehicle/ChSpindle.h"

#include "chrono_models/vehicle/generic/Generic_Wheel.h"
#include "chrono_models/vehicle/generic/tire/Generic_RigidTire.h"
#include "chrono_models/vehicle/marder/Marder.h"

using namespace chrono;
using namespace chrono::vehicle;

// -----------------------------------------------------------------------------

// Number of triangle-mesh visual shapes (with the given name, if not empty) in the visual model of a body.
static unsigned int CountMeshShapes(std::shared_ptr<ChBody> body, const std::string& name = "") {
    unsigned int count = 0;
    if (!body->GetVisualModel())
        return count;
    for (const auto& inst : body->GetVisualModel()->GetShapeInstances()) {
        auto mesh = std::dynamic_pointer_cast<ChVisualShapeTriangleMesh>(inst.shape);
        if (mesh && (name.empty() || mesh->GetName() == name))
            count++;
    }
    return count;
}

static std::string MissingMeshFile() {
    auto filename = (std::filesystem::temp_directory_path() / "chrono_utest_missing_mesh.obj").string();
    EXPECT_FALSE(std::filesystem::exists(filename));
    return filename;
}

static utils::ChBodyGeometry CreatePrimitiveGeometry() {
    utils::ChBodyGeometry geometry;
    geometry.vis_boxes.push_back(utils::ChBodyGeometry::BoxShape(ChVector3d(0, 0, 0), QUNIT, ChVector3d(1, 1, 1)));
    geometry.vis_spheres.push_back(utils::ChBodyGeometry::SphereShape(ChVector3d(1, 0, 0), 0.5));
    return geometry;
}

// -----------------------------------------------------------------------------

// A model file that does not exist: the primitive shapes are used instead (2 shapes: 1 box and 1 sphere).
TEST(BodyGeometryMeshVisualization, missing_model_file) {
    auto body = chrono_types::make_shared<ChBody>();

    auto geometry = CreatePrimitiveGeometry();
    geometry.vis_model_file = MissingMeshFile();

    geometry.CreateVisualizationAssets(body, VisualizationType::MESH);

    ASSERT_TRUE(body->GetVisualModel());
    EXPECT_EQ(body->GetVisualModel()->GetNumShapes(), 2u);
    EXPECT_EQ(CountMeshShapes(body), 0u);
}

// A model file that can be loaded: the mesh is the only visual shape (primitives are not used).
TEST(BodyGeometryMeshVisualization, valid_model_file) {
    // Single triangle
    auto filename = (std::filesystem::temp_directory_path() / "chrono_utest_triangle.obj").string();
    {
        std::ofstream obj(filename);
        obj << "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n";
    }

    auto body = chrono_types::make_shared<ChBody>();
    auto geometry = CreatePrimitiveGeometry();
    geometry.vis_model_file = filename;

    geometry.CreateVisualizationAssets(body, VisualizationType::MESH);
    std::filesystem::remove(filename);

    ASSERT_TRUE(body->GetVisualModel());
    EXPECT_EQ(body->GetVisualModel()->GetNumShapes(), 1u);
    EXPECT_EQ(CountMeshShapes(body, "chrono_utest_triangle"), 1u);
}

// A mesh shape created from a file that does not exist is skipped (no visual shape, no crash).
TEST(BodyGeometryMeshVisualization, missing_mesh_shape) {
    auto missing = MissingMeshFile();

    auto body = chrono_types::make_shared<ChBody>();
    utils::ChBodyGeometry geometry;
    geometry.vis_meshes.push_back(utils::ChBodyGeometry::TrimeshShape(VNULL, QUNIT, missing, 1.0));
    EXPECT_FALSE(geometry.vis_meshes[0].trimesh);

    geometry.CreateVisualizationAssets(body, VisualizationType::PRIMITIVES);
    EXPECT_TRUE(!body->GetVisualModel() || body->GetVisualModel()->GetNumShapes() == 0u);
}

// -----------------------------------------------------------------------------

// A null mesh is ignored: the shape keeps its (empty) mesh, or the mesh set before.
TEST(VisualShapeTriangleMesh, null_mesh) {
    ChVisualShapeTriangleMesh shape;
    ASSERT_TRUE(shape.GetMesh());

    EXPECT_NO_THROW(shape.SetMesh(nullptr));
    ASSERT_TRUE(shape.GetMesh());
    EXPECT_EQ(shape.GetMesh()->GetNumTriangles(), 0);

    auto mesh = chrono_types::make_shared<ChTriangleMeshConnected>();
    mesh->m_vertices = {ChVector3d(0, 0, 0), ChVector3d(1, 0, 0), ChVector3d(0, 1, 0)};
    mesh->m_face_v_indices = {ChVector3i(0, 1, 2)};
    shape.SetMesh(mesh);
    shape.SetMesh(nullptr);
    EXPECT_EQ(shape.GetMesh(), mesh);
    EXPECT_EQ(shape.GetMesh()->GetNumTriangles(), 1);
}

// A collision mesh that could not be loaded is skipped (collision shapes, bounding box, combination of geometries).
TEST(BodyGeometryMeshCollision, missing_mesh_shape) {
    auto body = chrono_types::make_shared<ChBody>();
    ChSystemNSC sys;
    sys.AddBody(body);

    utils::ChBodyGeometry geometry;
    geometry.materials.push_back(ChContactMaterialData());
    geometry.coll_boxes.push_back(utils::ChBodyGeometry::BoxShape(ChVector3d(0, 0, 0), QUNIT, ChVector3d(1, 2, 3), 0));
    geometry.coll_meshes.push_back(utils::ChBodyGeometry::TrimeshShape(VNULL, QUNIT, MissingMeshFile(), 1.0, 0.0, 0));
    ASSERT_FALSE(geometry.coll_meshes[0].trimesh);
    EXPECT_TRUE(geometry.HasCollision());

    // The bounding box is that of the box alone: [-0.5, 0.5] x [-1, 1] x [-1.5, 1.5]
    auto aabb = geometry.CalculateAABB();
    EXPECT_NEAR(aabb.min.x(), -0.5, 1e-12);
    EXPECT_NEAR(aabb.max.x(), +0.5, 1e-12);
    EXPECT_NEAR(aabb.min.y(), -1.0, 1e-12);
    EXPECT_NEAR(aabb.max.y(), +1.0, 1e-12);
    EXPECT_NEAR(aabb.min.z(), -1.5, 1e-12);
    EXPECT_NEAR(aabb.max.z(), +1.5, 1e-12);

    // Visualization of the collision shapes and collision shape creation skip the mesh
    EXPECT_NO_THROW(geometry.CreateVisualizationAssets(body, VisualizationType::COLLISION));
    ASSERT_TRUE(body->GetVisualModel());
    EXPECT_EQ(body->GetVisualModel()->GetNumShapes(), 1u);
    EXPECT_NO_THROW(geometry.CreateCollisionShapes(body, 1, ChContactMethod::NSC));

    // Combining geometries drops the mesh
    auto combined = utils::ChBodyGeometry::Combine({{ChFramed(), geometry}});
    EXPECT_EQ(combined->coll_boxes.size(), 1u);
    EXPECT_TRUE(combined->coll_meshes.empty());
}

// A terrain mesh file that cannot be loaded throws an exception that names the file; no patch is added.
TEST(RigidTerrainMesh, missing_mesh_file) {
    ChSystemNSC sys;
    RigidTerrain terrain(&sys);
    auto mat = ChContactMaterialData().CreateMaterial(ChContactMethod::NSC);

    auto filename = MissingMeshFile();
    try {
        terrain.AddPatch(mat, ChCoordsys<>(), filename, true, 0.01, false);
        FAIL() << "AddPatch did not throw";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find(filename), std::string::npos);
    }
    EXPECT_TRUE(terrain.GetPatches().empty());
}

TEST(SCMTerrainMesh, missing_mesh_file) {
    ChSystemNSC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    SCMTerrain terrain(&sys);
    EXPECT_THROW(terrain.Initialize(MissingMeshFile(), 0.1), std::invalid_argument);
}

// Rigid tire with access to the (protected) initialization function
class TestRigidTire : public generic::Generic_RigidTire {
  public:
    TestRigidTire(const std::string& name) : generic::Generic_RigidTire(name) {}
    void InitializeTire(std::shared_ptr<ChWheel> wheel) { Initialize(wheel); }
};

// A tire contact mesh file that cannot be loaded throws an exception that names the file.
TEST(RigidTireMesh, missing_contact_mesh_file) {
    ChSystemNSC sys;

    auto spindle = chrono_types::make_shared<ChSpindle>();
    sys.AddBody(spindle);
    auto wheel = chrono_types::make_shared<generic::Generic_Wheel>("wheel");
    wheel->Initialize(nullptr, spindle, VehicleSide::LEFT);

    auto filename = MissingMeshFile();
    TestRigidTire tire("tire");
    tire.SetContactMesh(filename, 0.005);
    try {
        tire.InitializeTire(wheel);
        FAIL() << "Initialize did not throw";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find(filename), std::string::npos);
    }
}

// -----------------------------------------------------------------------------

// Create and initialize the Marder (single-pin shoes) with mesh visualization.
// The track shoe, sprocket, idler, wheel and chassis meshes must all be found in the data directory.
TEST(MarderMeshVisualization, initialize) {
    ChSystemNSC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, -9.81));

    marder::Marder marder(&sys);
    marder.SetContactMethod(ChContactMethod::NSC);
    marder.SetChassisFixed(true);
    marder.SetInitPosition(ChCoordsys<>(ChVector3d(0, 0, 1), QUNIT));
    marder.Initialize();

    marder.SetChassisVisualizationType(VisualizationType::MESH);
    marder.SetSprocketVisualizationType(VisualizationType::MESH);
    marder.SetIdlerVisualizationType(VisualizationType::MESH);
    marder.SetIdlerWheelVisualizationType(VisualizationType::MESH);
    marder.SetRoadWheelVisualizationType(VisualizationType::MESH);
    marder.SetRollerVisualizationType(VisualizationType::MESH);
    marder.SetSuspensionVisualizationType(VisualizationType::PRIMITIVES);
    marder.SetTrackShoeVisualizationType(VisualizationType::MESH);

    auto& vehicle = marder.GetVehicle();
    for (auto side : {VehicleSide::LEFT, VehicleSide::RIGHT}) {
        auto track = vehicle.GetTrackAssembly(side);
        ASSERT_GT(track->GetNumTrackShoes(), 0u);

        // Every track shoe carries the shoe mesh (a loaded mesh replaces the primitive shapes)
        for (size_t i = 0; i < track->GetNumTrackShoes(); i++)
            EXPECT_EQ(CountMeshShapes(track->GetTrackShoe(i)->GetShoeBody(), "TrackShoeSinglePin"), 1u);
    }
    EXPECT_EQ(CountMeshShapes(marder.GetChassisBody(), "Hull"), 1u);

    // Assemble and advance the system
    sys.DoAssembly(AssemblyAnalysis::Level::FULL);
    EXPECT_NO_THROW(sys.DoStepDynamics(1e-3));
}
