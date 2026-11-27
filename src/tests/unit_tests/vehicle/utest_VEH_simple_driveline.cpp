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
// Tests for the kinematics and the torque transfer of the simple drivelines and of the 6WD shafts driveline:
// - ChSimpleDriveline: driveshaft speed = (alpha * w_front + (1 - alpha) * w_rear) / conical ratio, positive for
//   forward motion, i.e. the same as for a shafts-based driveline with the same conical gear ratio
// - ChSimpleDriveline: driveshaft torque x driveshaft speed = sum of wheel torque x wheel speed (no limited slip)
// - ChSimpleDriveline: the driveshaft speed is updated also when the driveshaft torque is within the zero-torque band,
//   and such a torque is not transmitted
// - ChSimpleDrivelineXWD: driveshaft speed = mean wheel speed / conical ratio, as for the 8WD shafts driveline
// - ChShaftsDriveline6WD: the front conical gear uses the front ratio and both rear conical gears the rear ratio
//
// =============================================================================

#include <cmath>
#include <memory>
#include <vector>

#include "gtest/gtest.h"

#include "chrono_vehicle/terrain/RigidTerrain.h"
#include "chrono_vehicle/wheeled_vehicle/driveline/ChSimpleDriveline.h"

#include "chrono_models/vehicle/hmmwv/HMMWV.h"
#include "chrono_models/vehicle/mrole/mrole.h"
#include "chrono_models/vehicle/mrole/mrole_VehicleFull.h"
#include "chrono_models/vehicle/mrole/driveline/mrole_Driveline6WD.h"

using namespace chrono;
using namespace chrono::vehicle;

static const double v0 = 10.0;    // initial forward speed (m/s)
static const double step = 1e-3;  // integration step size (s)

// -----------------------------------------------------------------------------

// HMMWV with the given driveline and simple-map power train, at v0 with wheel angular velocities omega.
// HMMWV_SimpleDriveline and HMMWV_Driveline4WD both have conical gear ratios 0.2; the front torque fraction is 0.5.
static std::shared_ptr<hmmwv::HMMWV_Full> CreateHMMWV(DrivelineTypeWV drive_type, double omega) {
    auto hmmwv = chrono_types::make_shared<hmmwv::HMMWV_Full>();
    hmmwv->SetContactMethod(ChContactMethod::NSC);
    hmmwv->SetDriveType(drive_type);
    hmmwv->SetEngineType(EngineModelType::SIMPLE_MAP);
    hmmwv->SetTransmissionType(TransmissionModelType::AUTOMATIC_SIMPLE_MAP);
    hmmwv->SetTireType(TireModelType::TMEASY);
    hmmwv->SetTireStepSize(step);
    hmmwv->SetInitPosition(ChCoordsys<>(ChVector3d(0, 0, 0.5), QUNIT));
    hmmwv->SetInitFwdVel(v0);
    hmmwv->SetInitWheelAngVel({omega, omega, omega, omega});
    hmmwv->Initialize();
    return hmmwv;
}

// Rigid terrain patch under the tires of the given vehicle.
static std::shared_ptr<RigidTerrain> CreateTerrain(ChWheeledVehicle& vehicle) {
    double ground = 1e9;
    for (auto& axle : vehicle.GetAxles())
        for (auto& wheel : axle->GetWheels())
            ground = std::min(ground, wheel->GetPos().z() - wheel->GetTire()->GetRadius());
    auto terrain = chrono_types::make_shared<RigidTerrain>(vehicle.GetSystem());
    ChContactMaterialData minfo;
    minfo.mu = 0.9f;
    auto mat = minfo.CreateMaterial(vehicle.GetSystem()->GetContactMethod());
    terrain->AddPatch(mat, ChCoordsys<>(ChVector3d(0, 0, ground + 0.015), QUNIT), 200, 20);
    terrain->Initialize();
    return terrain;
}

// Driveshaft speed of the HMMWV simple driveline from the axle speeds (axle shafts spin opposite to the wheels).
static double SimpleDriveshaftSpeed(ChWheeledVehicle& vehicle) {
    double front = -0.5 * (vehicle.GetSpindleOmega(0, LEFT) + vehicle.GetSpindleOmega(0, RIGHT));
    double rear = -0.5 * (vehicle.GetSpindleOmega(1, LEFT) + vehicle.GetSpindleOmega(1, RIGHT));
    return (0.5 * front + 0.5 * rear) / 0.2;
}

// -----------------------------------------------------------------------------

TEST(ChSimpleDriveline, SameSpeedAsShaftsDriveline) {
    double omega = v0 / 0.4699;
    auto simple = CreateHMMWV(DrivelineTypeWV::SIMPLE, omega);
    auto shafts = CreateHMMWV(DrivelineTypeWV::AWD, omega);

    double driveshaft = omega / 0.2;
    EXPECT_NEAR(simple->GetVehicle().GetDriveline()->GetOutputDriveshaftSpeed(), driveshaft, 1e-9 * driveshaft);
    EXPECT_NEAR(shafts->GetVehicle().GetDriveline()->GetOutputDriveshaftSpeed(), driveshaft, 1e-9 * driveshaft);

    auto pt_simple = simple->GetVehicle().GetPowertrainAssembly();
    auto pt_shafts = shafts->GetVehicle().GetPowertrainAssembly();
    EXPECT_EQ(pt_simple->GetTransmission()->GetCurrentGear(), pt_shafts->GetTransmission()->GetCurrentGear());
    EXPECT_NEAR(pt_simple->GetEngine()->GetMotorSpeed(), pt_shafts->GetEngine()->GetMotorSpeed(), 1e-9 * driveshaft);
}

TEST(ChSimpleDriveline, PowerBalance) {
    double omega = v0 / 0.4699;
    auto hmmwv = CreateHMMWV(DrivelineTypeWV::SIMPLE, omega);
    auto& vehicle = hmmwv->GetVehicle();
    auto terrain = CreateTerrain(vehicle);
    auto powertrain = vehicle.GetPowertrainAssembly();
    auto driveline = vehicle.GetDriveline();

    DriverInputs inputs = {0, 1, 0, 0};  // full throttle
    double time = 0;
    int checked = 0;
    for (int i = 0; i < 300; i++) {
        double torque = powertrain->GetOutputTorque();  // torque passed to the driveline in this Synchronize
        terrain->Synchronize(time);
        vehicle.Synchronize(time, inputs, *terrain);

        // Power balance holds if no differential limits the slip (left/right speed difference below 0.25 rad/s)
        bool open = true;
        double power_wheels = 0;
        for (int axle = 0; axle < 2; axle++) {
            open = open && std::abs(vehicle.GetSpindleOmega(axle, LEFT) - vehicle.GetSpindleOmega(axle, RIGHT)) < 0.25;
            for (auto side : {LEFT, RIGHT})
                power_wheels += driveline->GetSpindleTorque(axle, side) * (-vehicle.GetSpindleOmega(axle, side));
        }
        double power_driveshaft = (std::abs(torque) < 20 ? 0 : torque) * driveline->GetOutputDriveshaftSpeed();
        EXPECT_NEAR(driveline->GetOutputDriveshaftSpeed(), SimpleDriveshaftSpeed(vehicle), 1e-9);
        if (open && std::abs(power_driveshaft) > 0) {
            EXPECT_NEAR(power_wheels, power_driveshaft, 1e-9 * std::abs(power_driveshaft));
            checked++;
        }

        terrain->Advance(step);
        vehicle.Advance(step);
        time += step;
    }
    EXPECT_GT(checked, 100);
    // Traction: the vehicle accelerates
    EXPECT_GT(vehicle.GetSpeed(), v0);
}

TEST(ChSimpleDriveline, ZeroTorqueBand) {
    double omega = v0 / 0.4699;
    auto hmmwv = CreateHMMWV(DrivelineTypeWV::SIMPLE, omega);
    auto& vehicle = hmmwv->GetVehicle();
    auto terrain = CreateTerrain(vehicle);
    auto driveline = std::dynamic_pointer_cast<ChSimpleDriveline>(vehicle.GetDriveline());
    ASSERT_TRUE(driveline);
    driveline->SetZeroTorqueBand(1e9);  // no driveshaft torque is transmitted

    DriverInputs inputs = {0, 1, 0, 0};
    double time = 0;
    for (int i = 0; i < 200; i++) {
        terrain->Synchronize(time);
        vehicle.Synchronize(time, inputs, *terrain);
        EXPECT_NEAR(driveline->GetOutputDriveshaftSpeed(), SimpleDriveshaftSpeed(vehicle), 1e-9);
        for (int axle = 0; axle < 2; axle++)
            for (auto side : {LEFT, RIGHT})
                EXPECT_EQ(driveline->GetSpindleTorque(axle, side), 0);
        terrain->Advance(step);
        vehicle.Advance(step);
        time += step;
    }
    // No traction: the vehicle coasts
    EXPECT_LT(vehicle.GetSpeed(), v0);
}

// -----------------------------------------------------------------------------

TEST(ChSimpleDrivelineXWD, SameSpeedAsShaftsDriveline) {
    // mrole_SimpleDrivelineXWD and mrole_Driveline8WD: conical gear ratio 0.1 on each of the 4 driven axles
    double omega = 15.0;
    double driveshaft = omega / 0.1;
    for (auto type : {DrivelineTypeWV::SIMPLE_XWD, DrivelineTypeWV::AWD8}) {
        mrole::mrole_Full mrole;
        mrole.SetContactMethod(ChContactMethod::NSC);
        mrole.SetDriveType(type);
        mrole.SetTireType(TireModelType::TMEASY);
        mrole.SetInitPosition(ChCoordsys<>(ChVector3d(0, 0, 0.8), QUNIT));
        mrole.SetInitFwdVel(v0);
        mrole.SetInitWheelAngVel({omega, omega, omega, omega});
        mrole.Initialize();
        EXPECT_NEAR(mrole.GetVehicle().GetDriveline()->GetOutputDriveshaftSpeed(), driveshaft, 1e-9 * driveshaft);
    }
}

// -----------------------------------------------------------------------------

// 6WD driveline with different front and rear conical gear ratios (front 0.25, rear 0.2).
class Test_Driveline6WD : public mrole::mrole_Driveline6WD {
  public:
    Test_Driveline6WD() : mrole::mrole_Driveline6WD("Driveline") {}
    virtual double GetFrontConicalGearRatio() const override { return 0.25; }
};

// mrole vehicle (driven axles 0, 2, 3) with the test 6WD driveline.
class Test_Vehicle6WD : public mrole::mrole_VehicleFull {
  public:
    Test_Vehicle6WD(bool fixed)
        : mrole::mrole_VehicleFull(fixed, DrivelineTypeWV::AWD6, BrakeType::SHAFTS, SteeringTypeWV::PITMAN_ARM, false, ChContactMethod::NSC, CollisionType::NONE) {
        m_driveline = chrono_types::make_shared<Test_Driveline6WD>();
    }
};

TEST(ChShaftsDriveline6WD, ConicalGearRatios) {
    // Kinematics: the central differential averages the two rear conical gear inputs, w_rear / 0.2 each.
    {
        Test_Vehicle6WD vehicle(false);
        double omega = 15.0;
        vehicle.SetInitWheelAngVel({omega, omega, omega, omega});
        vehicle.Initialize(ChCoordsys<>(ChVector3d(0, 0, 1), QUNIT), v0);
        double driveshaft = omega / 0.2;
        EXPECT_NEAR(vehicle.GetDriveline()->GetOutputDriveshaftSpeed(), driveshaft, 1e-9 * driveshaft);
    }

    // Static torque split (fixed chassis, wheels held by the brakes): the central differential passes half of the
    // driveshaft torque to each rear conical gear, which multiplies it by 1 / 0.2. The front axle is driven only
    // through the hang-on clutch, which is open.
    {
        Test_Vehicle6WD vehicle(true);
        vehicle.Initialize(ChCoordsys<>(ChVector3d(0, 0, 1), QUNIT), 0);
        auto driveline = vehicle.GetDriveline();
        const double torque = 100;
        DriverInputs inputs = {0, 0, 1, 0};  // full braking
        double time = 0;
        for (int i = 0; i < 500; i++) {
            vehicle.Synchronize(time, inputs);
            driveline->Synchronize(time, inputs, torque);
            vehicle.Advance(step);
            time += step;
        }
        auto axle_torque = [&](int axle) { return driveline->GetSpindleTorque(axle, LEFT) + driveline->GetSpindleTorque(axle, RIGHT); };
        EXPECT_NEAR(axle_torque(2), 0.5 * torque / 0.2, 1e-2 * torque);
        EXPECT_NEAR(axle_torque(3), 0.5 * torque / 0.2, 1e-2 * torque);
        EXPECT_NEAR(axle_torque(0), 0, 1e-2 * torque);
    }
}
