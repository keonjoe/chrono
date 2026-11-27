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
// Authors: Radu Serban
// =============================================================================
//
// Simple driveline model. This template can be used to model a FWD, a RWD, or a
// 4WD driveline. It uses a constant front/rear torque split (a value between 0
// and 1 for 4WD, 0 for RWD, and 1 for FWD) and a simple model for a Torsen
// limited-slip differential (front and/or rear).
//
// =============================================================================

#include <cmath>

#include "chrono_vehicle/wheeled_vehicle/driveline/ChSimpleDriveline.h"

namespace chrono {
namespace vehicle {

// -----------------------------------------------------------------------------
// Construct a default 4WD simple driveline.
// -----------------------------------------------------------------------------
ChSimpleDriveline::ChSimpleDriveline(const std::string& name)
    : ChDrivelineWV(name),
      m_connected(true),
      m_front_locked(false),
      m_rear_locked(false),
      m_central_locked(false),
      m_driveshaft_speed(0),
      m_zero_torque(20) {}

// Torque bias ratio used to model a locked differential.
static const double kLockedBias = 1e3;

void ChSimpleDriveline::LockAxleDifferential(int axle, bool lock) {
    if (axle == 0 || axle == -1)
        m_front_locked = lock;
    if (axle == 1 || axle == -1)
        m_rear_locked = lock;
}

void ChSimpleDriveline::LockCentralDifferential(int which, bool lock) {
    m_central_locked = lock;
}

// -----------------------------------------------------------------------------
// Initialize the driveline subsystem.
// This function connects this driveline to the specified axles.
// -----------------------------------------------------------------------------
void ChSimpleDriveline::Initialize(std::shared_ptr<ChChassis> chassis,
                                   const ChAxleList& axles,
                                   const std::vector<int>& driven_axles) {
    assert(driven_axles.size() == 2);

    m_driven_axles = driven_axles;

    // Grab handles to the suspension wheel shafts.
    m_front_left = axles[m_driven_axles[0]]->m_suspension->GetAxle(LEFT);
    m_front_right = axles[m_driven_axles[0]]->m_suspension->GetAxle(RIGHT);

    m_rear_left = axles[m_driven_axles[1]]->m_suspension->GetAxle(LEFT);
    m_rear_right = axles[m_driven_axles[1]]->m_suspension->GetAxle(RIGHT);

    // Initialize the driveshaft speed from the initial wheel angular velocities (as calculated in Synchronize).
    double alpha = GetFrontTorqueFraction();
    double speed_front = -0.5 * (m_front_left->GetPosDt() + m_front_right->GetPosDt()) / GetFrontConicalGearRatio();
    double speed_rear = -0.5 * (m_rear_left->GetPosDt() + m_rear_right->GetPosDt()) / GetRearConicalGearRatio();
    m_driveshaft_speed = alpha * speed_front + (1 - alpha) * speed_rear;

    ChPart::Initialize();
}

// -----------------------------------------------------------------------------
// This utility function implements a simple model of Torsen limited-slip
// differential with a max_bias:1 torque bias ratio.
// We hard-code the speed difference range over which the torque bias grows from
// a value of 1 to a value of max_bias to the interval [0.25, 0.5].
// -----------------------------------------------------------------------------
void differentialSplit(double torque,
                       double max_bias,
                       double speed_left,
                       double speed_right,
                       double& torque_left,
                       double& torque_right) {
    double diff = std::abs(speed_left - speed_right);

    // The bias grows from 1 at diff=0.25 to max_bias at diff=0.5
    double bias = 1;
    if (diff > 0.5)
        bias = max_bias;
    else if (diff > 0.25)
        bias = 4 * (max_bias - 1) * diff + (2 - max_bias);

    // Split torque to the slow and fast wheels.
    double alpha = bias / (1 + bias);
    double slow = alpha * torque;
    double fast = torque - slow;

    if (std::abs(speed_left) < std::abs(speed_right)) {
        torque_left = slow;
        torque_right = fast;
    } else {
        torque_left = fast;
        torque_right = slow;
    }
}

// -----------------------------------------------------------------------------

void ChSimpleDriveline::Synchronize(double time, const DriverInputs& driver_inputs, double driveshaft_torque) {
    if (!m_connected)
        return;

    double alpha = GetFrontTorqueFraction();

    // Set driveshaft speed (output to transmission).
    // Each conical gear reduces the driveshaft speed by its gear ratio (w_differential = ratio * w_driveshaft, with the
    // differential speed the mean of the two axle speeds) and amplifies the torque by the inverse ratio. The axle shafts
    // spin opposite to the wheels and are driven by the negative of the wheel torque, so that a positive driveshaft
    // speed and torque correspond to forward motion and traction (as in the shafts-based drivelines). With the
    // driveshaft speed weighted by the torque fractions, the driveshaft power equals the power delivered to the axles.
    double speed_front = -0.5 * (m_front_left->GetPosDt() + m_front_right->GetPosDt()) / GetFrontConicalGearRatio();
    double speed_rear = -0.5 * (m_rear_left->GetPosDt() + m_rear_right->GetPosDt()) / GetRearConicalGearRatio();

    m_driveshaft_speed = alpha * speed_front + (1 - alpha) * speed_rear;

    // Driveshaft torques within the zero-torque band are not transmitted
    if (std::abs(driveshaft_torque) < m_zero_torque)
        driveshaft_torque = 0;

    // Split the input torque front/back (fixed fraction, or biased towards the slower axle if the central differential
    // is locked).
    double torque_front = driveshaft_torque * alpha;
    double torque_rear = driveshaft_torque - torque_front;
    if (m_central_locked)
        differentialSplit(driveshaft_torque, kLockedBias, speed_front, speed_rear, torque_front, torque_rear);

    torque_front /= GetFrontConicalGearRatio();
    torque_rear /= GetRearConicalGearRatio();

    // Split the axle torques for the corresponding left/right wheels and apply
    // them to the suspension wheel shafts.
    double torque_left;
    double torque_right;

    differentialSplit(torque_front, m_front_locked ? kLockedBias : GetFrontDifferentialMaxBias(),
                      m_front_left->GetPosDt(), m_front_right->GetPosDt(), torque_left, torque_right);
    m_front_left->SetAppliedLoad(-torque_left);
    m_front_right->SetAppliedLoad(-torque_right);

    differentialSplit(torque_rear, m_rear_locked ? kLockedBias : GetRearDifferentialMaxBias(), m_rear_left->GetPosDt(),
                      m_rear_right->GetPosDt(), torque_left, torque_right);
    m_rear_left->SetAppliedLoad(-torque_left);
    m_rear_right->SetAppliedLoad(-torque_right);
}

// -----------------------------------------------------------------------------

double ChSimpleDriveline::GetSpindleTorque(int axle, VehicleSide side) const {
    if (!m_connected)
        return 0;

    if (axle == m_driven_axles[0]) {
        switch (side) {
            case LEFT:
                return -m_front_left->GetAppliedLoad();
            case RIGHT:
                return -m_front_right->GetAppliedLoad();
        }
    } else if (axle == m_driven_axles[1]) {
        switch (side) {
            case LEFT:
                return -m_rear_left->GetAppliedLoad();
            case RIGHT:
                return -m_rear_right->GetAppliedLoad();
        }
    }

    return 0;
}

// -----------------------------------------------------------------------------
void ChSimpleDriveline::Disconnect() {
    m_connected = false;
}

}  // end namespace vehicle
}  // end namespace chrono
