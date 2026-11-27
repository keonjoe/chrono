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
// Tests for wheeled vehicles initialized with a non-zero forward speed v0 and wheel angular velocities v0/R:
// - all vehicle bodies start with the chassis velocity
// - the driveline, transmission and engine shafts start at speeds consistent with the wheel spin
//   (driveshaft = w / conical ratio, motorshaft = driveshaft / gear ratio)
// - an automatic transmission in automatic shift mode starts in the gear selected by its shift schedule
// - with zero driver inputs on flat terrain, the vehicle momentum does not jump in the first step and the speed loss
//   over 0.5 s is only due to rolling resistance and engine drag
// - JSON-specified vehicles accept initial wheel angular velocities
//
// =============================================================================

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "chrono_vehicle/ChVehicleDataPath.h"
#include "chrono_vehicle/terrain/RigidTerrain.h"
#include "chrono_vehicle/utils/ChVehicleUtilsJSON.h"
#include "chrono_vehicle/wheeled_vehicle/vehicle/WheeledVehicle.h"

#include "chrono_models/vehicle/hmmwv/HMMWV.h"
#include "chrono_models/vehicle/man/MAN_5t.h"
#include "chrono_models/vehicle/sedan/Sedan.h"

using namespace chrono;
using namespace chrono::vehicle;

static const double v0 = 10.0;    // initial forward speed (m/s)
static const double step = 1e-3;  // integration step size (s)

// -----------------------------------------------------------------------------

// Velocity of the center of mass of all bodies of the given vehicle (global X component).
static double ComSpeed(ChWheeledVehicle& vehicle) {
    double mass = 0;
    double momentum = 0;
    for (auto& body : vehicle.GetSystem()->GetBodies()) {
        if (VehicleObjTag::ExtractVehicleTag(body->GetTag()) != vehicle.GetVehicleTag())
            continue;
        mass += body->GetMass();
        momentum += body->GetMass() * body->GetPosDt().x();
    }
    return momentum / mass;
}

// Check that all vehicle bodies move with the chassis.
static void CheckBodyVelocities(ChWheeledVehicle& vehicle) {
    for (auto& body : vehicle.GetSystem()->GetBodies()) {
        if (VehicleObjTag::ExtractVehicleTag(body->GetTag()) != vehicle.GetVehicleTag())
            continue;
        EXPECT_NEAR(body->GetPosDt().x(), v0, 1e-12) << body->GetName();
        EXPECT_NEAR(body->GetPosDt().y(), 0, 1e-12) << body->GetName();
        EXPECT_NEAR(body->GetPosDt().z(), 0, 1e-12) << body->GetName();
    }
}

// Check the drivetrain speeds at initialization: driveshaft = w / conical_ratio, motorshaft = driveshaft / gear_ratio,
// engine speed = motorshaft speed (clamped to max_engine_speed for engine models that limit their speed).
static void CheckDrivetrainSpeeds(ChWheeledVehicle& vehicle, double omega, double conical_ratio, double gear_ratio, double max_engine_speed = 1e30) {
    double driveshaft = omega / conical_ratio;
    double motorshaft = driveshaft / gear_ratio;
    auto powertrain = vehicle.GetPowertrainAssembly();
    EXPECT_NEAR(vehicle.GetDriveline()->GetOutputDriveshaftSpeed(), driveshaft, 1e-9 * driveshaft);
    EXPECT_NEAR(powertrain->GetTransmission()->GetOutputMotorshaftSpeed(), motorshaft, 1e-9 * motorshaft);
    EXPECT_NEAR(powertrain->GetEngine()->GetMotorSpeed(), std::min(motorshaft, max_engine_speed), 1e-9 * motorshaft);
}

// Place a rigid terrain patch under the tires (with a small static deflection), run 0.5 s with zero driver inputs and
// check the vehicle momentum: no jump in the first step, and a speed loss explained by rolling resistance and engine
// drag (at most max_loss).
static void CheckCoasting(ChWheeledVehicle& vehicle, double max_loss) {
    double ground = 1e9;
    for (auto& axle : vehicle.GetAxles())
        for (auto& wheel : axle->GetWheels())
            ground = std::min(ground, wheel->GetPos().z() - wheel->GetTire()->GetRadius());

    RigidTerrain terrain(vehicle.GetSystem());
    ChContactMaterialData minfo;
    minfo.mu = 0.9f;
    minfo.cr = 0.01f;
    minfo.Y = 2e7f;
    auto mat = minfo.CreateMaterial(vehicle.GetSystem()->GetContactMethod());
    terrain.AddPatch(mat, ChCoordsys<>(ChVector3d(0, 0, ground + 0.015), QUNIT), 200, 20);
    terrain.Initialize();

    DriverInputs inputs = {0, 0, 0, 0};
    double time = 0;
    double v_com = ComSpeed(vehicle);
    EXPECT_NEAR(v_com, v0, 1e-12);

    for (int i = 0; i < 500; i++) {
        terrain.Synchronize(time);
        vehicle.Synchronize(time, inputs, terrain);
        terrain.Advance(step);
        vehicle.Advance(step);
        time += step;
        double v_new = ComSpeed(vehicle);
        if (i == 0) {
            // Before the drivetrain and the vehicle bodies were initialized consistently, the first step removed
            // 2 to 3 m/s (an acceleration of the order of -2000 m/s^2).
            EXPECT_LT(std::abs(v_new - v_com) / step, 0.5);
        }
        v_com = v_new;
    }

    EXPECT_LT(v_com, v0);
    EXPECT_GT(v_com, v0 - max_loss);
}

// -----------------------------------------------------------------------------

// HMMWV_Full: driveline, transmission and engine type; conical gear ratio 0.2 (HMMWV_Driveline2WD/4WD);
// HMMWV_AutomaticTransmissionShafts: gear ratios 0.2, 0.4, 0.8, upshift at 2500 rpm (261.8 rad/s). At v0 = 10 m/s the
// driveshaft speed is 106.4 rad/s, so the shift schedule selects 3rd gear (in-gear shaft speed 133.0 rad/s).
static void TestHMMWV(DrivelineTypeWV drive_type, double max_loss) {
    hmmwv::HMMWV_Full hmmwv;
    hmmwv.SetContactMethod(ChContactMethod::NSC);
    hmmwv.SetDriveType(drive_type);
    hmmwv.SetEngineType(EngineModelType::SHAFTS);
    hmmwv.SetTransmissionType(TransmissionModelType::AUTOMATIC_SHAFTS);
    hmmwv.SetTireType(TireModelType::TMEASY);
    hmmwv.SetTireStepSize(step);
    hmmwv.SetInitPosition(ChCoordsys<>(ChVector3d(0, 0, 0.5), QUNIT));
    hmmwv.SetInitFwdVel(v0);
    const double radius = 0.4699;  // HMMWV_TMeasyTire unloaded radius
    double omega = v0 / radius;
    hmmwv.SetInitWheelAngVel({omega, omega, omega, omega});
    hmmwv.Initialize();

    auto& vehicle = hmmwv.GetVehicle();
    ASSERT_NEAR(vehicle.GetTire(0, LEFT)->GetRadius(), radius, 1e-4);
    CheckBodyVelocities(vehicle);
    EXPECT_EQ(vehicle.GetPowertrainAssembly()->GetTransmission()->GetCurrentGear(), 3);
    CheckDrivetrainSpeeds(vehicle, omega, 0.2, 0.8);
    CheckCoasting(vehicle, max_loss);
}

TEST(ChVehicleInitialSpeed, HMMWV_4WD_shafts) {
    // Expected loss over 0.5 s: rolling resistance (0.01 g) and engine drag in 3rd gear (0.17 m/s in total)
    TestHMMWV(DrivelineTypeWV::AWD, 0.3);
}

TEST(ChVehicleInitialSpeed, HMMWV_2WD_shafts) {
    TestHMMWV(DrivelineTypeWV::RWD, 0.3);
}

TEST(ChVehicleInitialSpeed, Sedan) {
    // Sedan: ChShaftsDriveline2WD with conical gear ratio 0.2, simple-map engine and transmission.
    sedan::Sedan sedan;
    sedan.SetContactMethod(ChContactMethod::NSC);
    sedan.SetTireType(TireModelType::TMEASY);
    sedan.SetTireStepSize(step);
    sedan.SetInitPosition(ChCoordsys<>(ChVector3d(0, 0, 0.5), QUNIT));
    sedan.SetInitFwdVel(v0);
    const double radius = 0.3266;  // Sedan_TMeasyTire unloaded radius
    double omega = v0 / radius;
    sedan.SetInitWheelAngVel({omega, omega, omega, omega});
    sedan.Initialize();

    auto& vehicle = sedan.GetVehicle();
    ASSERT_NEAR(vehicle.GetTire(0, LEFT)->GetRadius(), radius, 1e-3);
    CheckBodyVelocities(vehicle);
    double driveshaft = omega / 0.2;
    EXPECT_NEAR(vehicle.GetDriveline()->GetOutputDriveshaftSpeed(), driveshaft, 1e-9 * driveshaft);
    // Simple-map transmission: motorshaft speed = |driveshaft speed| / gear ratio (first gear ratio < 1);
    // simple-map engine: engine speed = motorshaft speed, clamped to the maximum engine speed.
    auto powertrain = vehicle.GetPowertrainAssembly();
    EXPECT_GT(powertrain->GetTransmission()->GetOutputMotorshaftSpeed(), driveshaft);
    EXPECT_GT(powertrain->GetEngine()->GetMotorSpeed(), driveshaft);
    CheckCoasting(vehicle, 0.3);
}

TEST(ChVehicleInitialSpeed, MAN_5t) {
    // MAN_5t: ChShaftsDriveline4WD with conical gear ratio 1/(3.947*1.706*1.02), simple-map powertrain.
    man::MAN_5t truck;
    truck.SetContactMethod(ChContactMethod::NSC);
    truck.SetTireType(TireModelType::TMEASY);
    truck.SetTireStepSize(step);
    truck.SetInitPosition(ChCoordsys<>(ChVector3d(0, 0, 0.8), QUNIT));
    truck.SetInitFwdVel(v0);
    const double radius = 0.629;  // MAN_5t_TMeasyTire unloaded radius
    double omega = v0 / radius;
    truck.SetInitWheelAngVel({omega, omega, omega, omega});
    truck.Initialize();

    auto& vehicle = truck.GetVehicle();
    ASSERT_NEAR(vehicle.GetTire(0, LEFT)->GetRadius(), radius, 1e-3);
    CheckBodyVelocities(vehicle);
    double driveshaft = omega * 3.947 * 1.706 * 1.02;
    EXPECT_NEAR(vehicle.GetDriveline()->GetOutputDriveshaftSpeed(), driveshaft, 1e-9 * driveshaft);
    auto powertrain = vehicle.GetPowertrainAssembly();
    EXPECT_GT(powertrain->GetTransmission()->GetOutputMotorshaftSpeed(), driveshaft);
    CheckCoasting(vehicle, 0.3);
}

TEST(ChVehicleInitialSpeed, JSON_vehicle) {
    WheeledVehicle vehicle(GetVehicleDataFile("hmmwv/vehicle/HMMWV_Vehicle_4WD.json"), ChContactMethod::NSC, false, false);
    EXPECT_THROW(vehicle.SetInitWheelAngVel({1.0}), std::invalid_argument);

    const double omega = 20;
    vehicle.SetInitWheelAngVel(std::vector<double>(2 * vehicle.GetNumberAxles(), omega));
    vehicle.Initialize(ChCoordsys<>(ChVector3d(0, 0, 0.5), QUNIT), v0);

    // Axle shafts spin opposite to the wheels (see ChSuspension::GetAxleSpeed)
    for (unsigned int i = 0; i < vehicle.GetNumberAxles(); i++) {
        EXPECT_NEAR(vehicle.GetSpindleOmega(i, LEFT), -omega, 1e-12);
        EXPECT_NEAR(vehicle.GetSpindleOmega(i, RIGHT), -omega, 1e-12);
    }
    CheckBodyVelocities(vehicle);

    // Powertrain created after the vehicle: shafts engine and automatic shafts transmission (gear ratios 0.2, 0.4, 0.8,
    // upshift at 261.8 rad/s), driveline HMMWV_Driveline4WD.json (conical gear ratio 0.2). The driveshaft speed is
    // 100 rad/s, so the shift schedule selects 2nd gear.
    auto engine = ReadEngineJSON(GetVehicleDataFile("hmmwv/powertrain/HMMWV_EngineShafts.json"));
    auto transmission = ReadTransmissionJSON(GetVehicleDataFile("hmmwv/powertrain/HMMWV_AutomaticTransmissionShafts.json"));
    vehicle.InitializePowertrain(chrono_types::make_shared<ChPowertrainAssembly>(engine, transmission));
    EXPECT_EQ(transmission->GetCurrentGear(), 2);
    CheckDrivetrainSpeeds(vehicle, omega, 0.2, 0.4);
}

TEST(ChVehicleInitialSpeed, AtRest) {
    // A vehicle initialized at rest keeps all drivetrain speeds at zero.
    hmmwv::HMMWV_Full hmmwv;
    hmmwv.SetEngineType(EngineModelType::SHAFTS);
    hmmwv.SetTransmissionType(TransmissionModelType::AUTOMATIC_SHAFTS);
    hmmwv.SetInitPosition(ChCoordsys<>(ChVector3d(0, 0, 0.5), QUNIT));
    hmmwv.Initialize();
    auto& vehicle = hmmwv.GetVehicle();
    EXPECT_EQ(vehicle.GetDriveline()->GetOutputDriveshaftSpeed(), 0);
    EXPECT_EQ(vehicle.GetPowertrainAssembly()->GetTransmission()->GetOutputMotorshaftSpeed(), 0);
    EXPECT_EQ(vehicle.GetPowertrainAssembly()->GetEngine()->GetMotorSpeed(), 0);
    for (auto& body : vehicle.GetSystem()->GetBodies())
        EXPECT_EQ(body->GetPosDt().Length(), 0) << body->GetName();
}
