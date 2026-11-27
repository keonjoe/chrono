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
// Unit test for Bullet ellipsoid collision: an ellipsoid, rotated about its X
// axis and pushed a prescribed depth into a box, must produce one contact with
// that penetration, a vertical normal and the analytic contact point, for all
// orientations (not only those with an ellipsoid axis along the normal).
//
// =============================================================================

#include <cmath>
#include <functional>
#include <string>

#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChContactContainer.h"

#include "gtest/gtest.h"

using namespace chrono;

class ContactReporter : public ChContactContainer::ReportContactCallback {
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
        normal = plane_coord.GetAxisX();
        // contact point on the ellipsoid (the lower of the two points)
        point = pA.z() < pB.z() ? pA : pB;
        return true;
    }

    int num_contacts = 0;
    double dist = 0;
    ChVector3d normal;
    ChVector3d point;
};

// Collide an ellipsoid with semi-axes (a, b, c), rotated by 'angle' about X, with a box whose top face is at z = 0.
static void CheckEllipsoidOnBox(double a, double b, double c, double angle, double envelope, double margin) {
    const double depth = 1e-4;

    // the system constructor sets the default envelope and margin: set them afterwards
    ChSystemSMC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    ChCollisionModel::SetDefaultSuggestedEnvelope(envelope);
    ChCollisionModel::SetDefaultSuggestedMargin(margin);
    auto mat = chrono_types::make_shared<ChContactMaterialSMC>();

    auto box = chrono_types::make_shared<ChBodyEasyBox>(2.0, 2.0, 0.2, 1000, false, true, mat);
    box->SetPos(ChVector3d(0, 0, -0.1));
    box->SetFixed(true);
    sys.Add(box);

    // height of the center above the lowest point and lateral position of that point
    double s = std::sin(angle);
    double co = std::cos(angle);
    double h = std::sqrt(b * b * s * s + c * c * co * co);
    double y_contact = s * co * (c * c - b * b) / h;

    auto ellipsoid = chrono_types::make_shared<ChBodyEasyEllipsoid>(ChVector3d(2 * a, 2 * b, 2 * c), 1000, false, true, mat);
    ellipsoid->SetPos(ChVector3d(0, 0, h - depth));
    ellipsoid->SetRot(QuatFromAngleX(angle));
    sys.Add(ellipsoid);

    sys.GetCollisionSystem()->BindAll();
    sys.ComputeCollisions();

    auto reporter = chrono_types::make_shared<ContactReporter>();
    sys.GetContactContainer()->ReportAllContacts(reporter);

    std::cout << "semi-axes (" << a << ", " << b << ", " << c << ")  angle " << angle * CH_RAD_TO_DEG << "  contacts " << reporter->num_contacts << "  distance " << reporter->dist
              << "  y " << reporter->point.y() << " (exact " << y_contact << ")" << std::endl;

    ASSERT_EQ(reporter->num_contacts, 1);
    ASSERT_NEAR(reporter->dist, -depth, 0.4 * depth);
    ASSERT_NEAR(std::abs(reporter->normal.z()), 1.0, 1e-4);
    ASSERT_NEAR(reporter->point.y(), y_contact, 0.03 * std::max(a, std::max(b, c)));
}

TEST(BulletCollision, EllipsoidOnBox) {
    const double envelopes[] = {0.03, 0.002};
    const double margins[] = {0.01, 0.001};
    for (int k = 0; k < 2; k++) {
        for (double deg = 0; deg <= 90; deg += 7.5) {
            // prolate, oblate, small and large
            CheckEllipsoidOnBox(0.03, 0.03, 0.06, deg * CH_DEG_TO_RAD, envelopes[k], margins[k]);
            CheckEllipsoidOnBox(0.06, 0.06, 0.03, deg * CH_DEG_TO_RAD, envelopes[k], margins[k]);
            CheckEllipsoidOnBox(0.2, 0.3, 0.6, deg * CH_DEG_TO_RAD, envelopes[k], margins[k]);
        }
    }
}

// Collide a fixed body A at the origin with an ellipsoid B (semi-axes sB, rotation qB) whose center is pushed a
// prescribed depth below the touching height zB along Z; the contact normal is along Z by symmetry.
static void CheckPairAlongZ(std::function<std::shared_ptr<ChBody>(std::shared_ptr<ChContactMaterial>)> makeA,
                            const ChVector3d& sB,
                            const ChQuaterniond& qB,
                            double zB,
                            double envelope,
                            double margin,
                            const std::string& name) {
    const double depth = 1e-3 * zB;

    // the system constructor sets the default envelope and margin: set them afterwards, before creating the bodies
    ChSystemSMC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    ChCollisionModel::SetDefaultSuggestedEnvelope(envelope);
    ChCollisionModel::SetDefaultSuggestedMargin(margin);
    auto mat = chrono_types::make_shared<ChContactMaterialSMC>();

    auto A = makeA(mat);
    A->SetFixed(true);
    sys.Add(A);

    auto B = chrono_types::make_shared<ChBodyEasyEllipsoid>(2.0 * sB, 1000, false, true, mat);
    B->SetPos(ChVector3d(0, 0, zB - depth));
    B->SetRot(qB);
    sys.Add(B);

    sys.GetCollisionSystem()->BindAll();
    sys.ComputeCollisions();

    auto reporter = chrono_types::make_shared<ContactReporter>();
    sys.GetContactContainer()->ReportAllContacts(reporter);

    std::cout << name << "  envelope " << envelope << "  contacts " << reporter->num_contacts << "  distance " << reporter->dist << " (exact " << -depth << ")  normal z "
              << reporter->normal.z() << std::endl;

    ASSERT_EQ(reporter->num_contacts, 1);
    // single precision GJK on two curved shapes is accurate to a few 1e-4 of the size (errors before the fix of the
    // degenerate tetrahedron test were 10 to 500 times the depth)
    ASSERT_NEAR(reporter->dist, -depth, 0.2 * depth);
    ASSERT_NEAR(std::abs(reporter->normal.z()), 1.0, 1e-3);
}

// Height of the center of an ellipsoid with semi-axes s and rotation q above its lowest point.
static double SupportHeight(const ChVector3d& s, const ChQuaterniond& q) {
    ChVector3d n = q.RotateBack(ChVector3d(0, 0, 1));
    return ChVector3d(s.x() * n.x(), s.y() * n.y(), s.z() * n.z()).Length();
}

TEST(BulletCollision, EllipsoidPairs) {
    using Mat = std::shared_ptr<ChContactMaterial>;
    const double envelopes[] = {0.03, 0.002};
    const double margins[] = {0.01, 0.001};
    for (int k = 0; k < 2; k++) {
        for (double scale : {1.0, 3.0, 10.0}) {
            const std::string sc = ", scale " + std::to_string(scale);
            // flat ellipsoid B, yawed about Z, stacked on a flat ellipsoid A
            for (double yaw : {0.0, 30.0, 60.0}) {
                auto A = [=](Mat mat) { return chrono_types::make_shared<ChBodyEasyEllipsoid>(scale * ChVector3d(0.12, 0.12, 0.06), 1000, false, true, mat); };
                CheckPairAlongZ(A, scale * ChVector3d(0.05, 0.04, 0.02), QuatFromAngleZ(yaw * CH_DEG_TO_RAD), scale * 0.05, envelopes[k], margins[k],
                                "ellipsoid-ellipsoid stack, yaw " + std::to_string(yaw) + sc);
            }
            // the same, on a sphere made from an ellipsoid shape and on a sphere shape
            auto S1 = [=](Mat mat) { return chrono_types::make_shared<ChBodyEasyEllipsoid>(scale * ChVector3d(0.2, 0.2, 0.2), 1000, false, true, mat); };
            CheckPairAlongZ(S1, scale * ChVector3d(0.05, 0.04, 0.02), QuatFromAngleZ(30 * CH_DEG_TO_RAD), scale * 0.12, envelopes[k], margins[k],
                            "ellipsoid-ellipsoid (round)" + sc);
            auto S2 = [=](Mat mat) { return chrono_types::make_shared<ChBodyEasySphere>(scale * 0.1, 1000, false, true, mat); };
            CheckPairAlongZ(S2, scale * ChVector3d(0.05, 0.04, 0.02), QuatFromAngleZ(30 * CH_DEG_TO_RAD), scale * 0.12, envelopes[k], margins[k], "ellipsoid-sphere" + sc);
            // mirrored tilted pair: B = A rotated by 180 degrees about Z, touching on the plane z = h
            for (auto angles : {ChVector3d(20, 35, 50), ChVector3d(30, 30, 30)}) {
                ChVector3d s = scale * ChVector3d(0.02, 0.03, 0.05);
                ChQuaterniond qA = QuatFromAngleZ(angles.z() * CH_DEG_TO_RAD) * QuatFromAngleY(angles.y() * CH_DEG_TO_RAD) * QuatFromAngleX(angles.x() * CH_DEG_TO_RAD);
                ChQuaterniond qB = QuatFromAngleZ(CH_PI) * qA;
                auto A = [=](Mat mat) {
                    auto body = chrono_types::make_shared<ChBodyEasyEllipsoid>(2.0 * s, 1000, false, true, mat);
                    body->SetRot(qA);
                    return body;
                };
                CheckPairAlongZ(A, s, qB, 2 * SupportHeight(s, qA), envelopes[k], margins[k], "mirrored ellipsoid pair" + sc);
            }
            // flat ellipsoid lying across a cylinder with its axis along Y (the cylinder is kept short: with single
            // precision Bullet, GJK loses accuracy for shapes much longer than the contact region, e.g. a cylinder of
            // length 8 radii gives errors of 5 to 25 % of a penetration of 1e-3 of the size)
            auto C = [=](Mat mat) { return chrono_types::make_shared<ChBodyEasyCylinder>(ChAxis::Y, scale * 0.15, scale * 0.6, 1000, false, true, mat); };
            CheckPairAlongZ(C, scale * ChVector3d(0.06, 0.04, 0.03), QUNIT, scale * 0.18, envelopes[k], margins[k], "ellipsoid-cylinder" + sc);
        }
    }
}
