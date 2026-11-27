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
// Unit test for Bullet contacts of a sphere with long shapes: a sphere resting on
// the top of a cylinder, a box, a rounded cylinder or a capsule lying along X (length 1, 10, 100 m)
// must get its contact point directly below its center and a vertical normal,
// at any position along the shape and any rotation of the sphere, independent of
// the length of the shape. An ellipsoid and a cylinder lying across the cylinder
// must get a normal without component along the axis of the cylinder (and, for the
// crossed cylinders, without component along either axis).
//
// =============================================================================

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>

#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChContactContainer.h"
#include "chrono/collision/ChCollisionShapeCapsule.h"
#include "chrono/collision/ChCollisionShapeRoundedCylinder.h"

#include "gtest/gtest.h"

using namespace chrono;

class SphereContactReporter : public ChContactContainer::ReportContactCallback {
  public:
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
        num_contacts++;
        dist = distance;
        point_A = pA;
        point_B = pB;
        normal = plane_coord.GetAxisX();
        return true;
    }

    int num_contacts = 0;
    double dist = 0;
    ChVector3d point_A;
    ChVector3d point_B;
    ChVector3d normal;
};

// Sphere of radius rs resting with a given depth on the top (z = h) of a fixed shape lying along X; the sphere is moved
// along the shape and rotated about Y as if it rolled.
static void CheckSphereOnLongShape(std::function<std::shared_ptr<ChBody>(std::shared_ptr<ChContactMaterial>)> makeShape, double length, double h, const std::string& name) {
    const double rs = 0.05;
    const double depth = 1e-6;

    ChSystemSMC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    ChCollisionModel::SetDefaultSuggestedMargin(0.001);
    auto mat = chrono_types::make_shared<ChContactMaterialSMC>();

    auto shape = makeShape(mat);
    shape->SetFixed(true);
    sys.Add(shape);

    auto sphere = chrono_types::make_shared<ChBodyEasySphere>(rs, 1000, false, true, mat);
    sys.Add(sphere);
    sys.GetCollisionSystem()->BindAll();

    // rounding of the sphere position along the shape in single precision Bullet
    const double tol = 1e-12 + 1e-7 * length;
    double worst_x = 0, worst_n = 0, worst_d = 0;
    int missing = 0;
    const int n = 1000;
    for (int i = 0; i < n; i++) {
        double x = (-0.45 + 0.9 * i / (n - 1.0)) * length;
        sphere->SetPos(ChVector3d(x, 0, h + rs - depth));
        sphere->SetRot(QuatFromAngleY(x / rs));
        sys.ComputeCollisions();
        auto reporter = chrono_types::make_shared<SphereContactReporter>();
        sys.GetContactContainer()->ReportAllContacts(reporter);
        if (reporter->num_contacts != 1) {
            missing++;
            continue;
        }
        worst_x = std::max({worst_x, std::abs(reporter->point_A.x() - x), std::abs(reporter->point_B.x() - x)});
        worst_n = std::max(worst_n, std::abs(reporter->normal.x()));
        worst_d = std::max(worst_d, std::abs(reporter->dist + depth));
    }

    std::cout << name << "  length " << length << "  positions without one contact " << missing << "  contact point x error " << worst_x << "  normal x " << worst_n
              << "  distance error " << worst_d << std::endl;

    ASSERT_EQ(missing, 0);
    ASSERT_LT(worst_x, tol);
    ASSERT_LT(worst_n, tol / rs);
    ASSERT_LT(worst_d, 1e-7);
}

TEST(BulletCollision, SphereOnLongShapes) {
    using Mat = std::shared_ptr<ChContactMaterial>;
    const double r = 0.1;
    for (double length : {1.0, 10.0, 100.0}) {
        auto cylinder = [=](Mat mat) { return chrono_types::make_shared<ChBodyEasyCylinder>(ChAxis::X, r, length, 1000, false, true, mat); };
        CheckSphereOnLongShape(cylinder, length, r, "cylinder");

        auto box = [=](Mat mat) { return chrono_types::make_shared<ChBodyEasyBox>(length, 2 * r, 2 * r, 1000, false, true, mat); };
        CheckSphereOnLongShape(box, length, r, "box");

        auto along_x = [](std::shared_ptr<ChCollisionShape> shape) {
            auto body = chrono_types::make_shared<ChBody>();
            body->AddCollisionShape(shape, ChFrame<>(VNULL, QuatFromAngleY(CH_PI_2)));
            body->EnableCollision(true);
            return body;
        };

        auto rounded_cylinder = [=](Mat mat) { return along_x(chrono_types::make_shared<ChCollisionShapeRoundedCylinder>(mat, r, length, 0.02)); };
        CheckSphereOnLongShape(rounded_cylinder, length, r, "rounded cylinder");

        auto capsule = [=](Mat mat) { return along_x(chrono_types::make_shared<ChCollisionShapeCapsule>(mat, r, length)); };
        CheckSphereOnLongShape(capsule, length, r, "capsule");
    }
}

// Worst normal components over all contacts
class NormalReporter : public ChContactContainer::ReportContactCallback {
  public:
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
        num_contacts++;
        ChVector3d n = plane_coord.GetAxisX();
        normal_x = std::max(normal_x, std::abs(n.x()));
        normal_y = std::max(normal_y, std::abs(n.y()));
        min_dist = std::min(min_dist, distance);
        return true;
    }

    int num_contacts = 0;
    double normal_x = 0;
    double normal_y = 0;
    double min_dist = 1e30;
};

// A body resting with a given depth on top of a fixed cylinder of radius r along X, moved along the cylinder and rotated
// about Y as if it rolled; 'height' gives the height of the lowest point of the body below its center for a rotation.
static void CheckBodyOnLongCylinder(std::function<std::shared_ptr<ChBody>(std::shared_ptr<ChContactMaterial>)> makeBody,
                                    std::function<double(double)> height,
                                    double roll_radius,
                                    double length,
                                    double tol_axial,
                                    double tol_lateral,
                                    const std::string& name) {
    const double r = 0.1;
    const double depth = 1e-5;

    ChSystemSMC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    ChCollisionModel::SetDefaultSuggestedEnvelope(0.001);
    ChCollisionModel::SetDefaultSuggestedMargin(0.001);
    auto mat = chrono_types::make_shared<ChContactMaterialSMC>();

    auto cylinder = chrono_types::make_shared<ChBodyEasyCylinder>(ChAxis::X, r, length, 1000, false, true, mat);
    cylinder->SetFixed(true);
    sys.Add(cylinder);
    auto body = makeBody(mat);
    sys.Add(body);
    sys.GetCollisionSystem()->BindAll();

    double worst_x = 0, worst_y = 0, worst_d = 0;
    int missing = 0;
    const int n = 1000;
    for (int i = 0; i < n; i++) {
        double x = (-0.45 + 0.9 * i / (n - 1.0)) * length;
        double angle = x / roll_radius;
        body->SetPos(ChVector3d(x, 0, r + height(angle) - depth));
        body->SetRot(QuatFromAngleY(angle));
        sys.ComputeCollisions();
        auto reporter = chrono_types::make_shared<NormalReporter>();
        sys.GetContactContainer()->ReportAllContacts(reporter);
        if (reporter->num_contacts == 0) {
            std::cout << "  no contact at x = " << x << std::endl;
            missing++;
            continue;
        }
        worst_x = std::max(worst_x, reporter->normal_x);
        worst_y = std::max(worst_y, reporter->normal_y);
        worst_d = std::max(worst_d, std::abs(reporter->min_dist + depth));
    }

    std::cout << name << "  length " << length << "  positions without contact " << missing << "  normal x " << worst_x << "  normal y " << worst_y << "  distance error "
              << worst_d << std::endl;

    ASSERT_EQ(missing, 0);
    ASSERT_LT(worst_x, tol_axial);
    ASSERT_LT(worst_y, tol_lateral);
}

// The normals of curved shapes against the lateral side of a cylinder: the component along the axis vanishes by
// construction; the lateral component of the ellipsoid normal is the accuracy of the GJK witness point relative to the
// cylinder radius (a few 1e-5 at 10 m with these margins), not rounding.
TEST(BulletCollision, CurvedShapesOnLongCylinder) {
    const ChVector3d semi(0.05, 0.04, 0.03);
    for (double length : {1.0, 10.0}) {
        auto ellipsoid = [=](std::shared_ptr<ChContactMaterial> mat) { return chrono_types::make_shared<ChBodyEasyEllipsoid>(2.0 * semi, 1000, false, true, mat); };
        auto height = [=](double a) { return std::sqrt(std::pow(semi.x() * std::sin(a), 2) + std::pow(semi.z() * std::cos(a), 2)); };
        CheckBodyOnLongCylinder(ellipsoid, height, semi.y(), length, 1e-12, 5e-4, "ellipsoid");
    }
    for (double length : {1.0, 10.0, 100.0}) {
        const double rc = 0.05;
        auto crossed = [=](std::shared_ptr<ChContactMaterial> mat) { return chrono_types::make_shared<ChBodyEasyCylinder>(ChAxis::Y, rc, 0.2, 1000, false, true, mat); };
        CheckBodyOnLongCylinder(crossed, [=](double) { return rc; }, rc, length, 1e-12, 1e-12, "crossed cylinder");
    }
}
