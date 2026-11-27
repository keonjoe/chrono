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
// Authors: Radu Serban, Rainer Gericke
// =============================================================================
//
// Simple driveline model. This template can be used to model a XWD driveline.
// Number of axles can be 1 to X.
// It uses a constant torque split depending on the number of axles driven and a
// simple model for Torsen limited-slip differentials.
//
// =============================================================================

#include <algorithm>
#include <cmath>

#include "chrono_vehicle/wheeled_vehicle/driveline/ChSimpleDrivelineXWD.h"

namespace chrono {
namespace vehicle {

// -----------------------------------------------------------------------------
// Construct a default 4WD simple driveline.
// -----------------------------------------------------------------------------
ChSimpleDrivelineXWD::ChSimpleDrivelineXWD(const std::string& name)
    : ChDrivelineWV(name),
      m_connected(true),
      m_all_axles_locked(false),
      m_central_locked(false),
      m_driveshaft_speed(0) {}

// Torque bias ratio used to model a locked differential.
static const double kLockedBiasXWD = 1e3;

void ChSimpleDrivelineXWD::LockAxleDifferential(int axle, bool lock) {
    if (axle < 0) {
        m_all_axles_locked = lock;
        m_axle_locked.assign(m_axle_locked.size(), lock);
        return;
    }
    if (axle >= (int)m_axle_locked.size())
        m_axle_locked.resize(axle + 1, m_all_axles_locked);
    m_axle_locked[axle] = lock;
}

void ChSimpleDrivelineXWD::LockCentralDifferential(int which, bool lock) {
    m_central_locked = lock;
}

// -----------------------------------------------------------------------------
// Initialize the driveline subsystem.
// This function connects this driveline to the specified axles.
// -----------------------------------------------------------------------------
void ChSimpleDrivelineXWD::Initialize(std::shared_ptr<ChChassis> chassis,
                                      const ChAxleList& axles,
                                      const std::vector<int>& driven_axles) {
    assert(axles.size() >= 1);

    m_driven_axles = driven_axles;

    for (int i = 0; i < driven_axles.size(); i++) {
        // Grab handles to the suspension wheel shafts.
        m_shaft_left.push_back(axles[m_driven_axles[i]]->m_suspension->GetAxle(LEFT));
        m_shaft_right.push_back(axles[m_driven_axles[i]]->m_suspension->GetAxle(RIGHT));
    }

    // Initialize the driveshaft speed from the initial wheel angular velocities (as calculated in Synchronize).
    double alpha = 1.0 / double(m_shaft_left.size());
    m_driveshaft_speed = 0;
    for (int i = 0; i < m_shaft_left.size(); i++) {
        m_driveshaft_speed -= alpha * 0.5 * (m_shaft_left[i]->GetPosDt() + m_shaft_right[i]->GetPosDt());
    }
    m_driveshaft_speed /= GetConicalGearRatio();

    ChPart::Initialize();
}

// -----------------------------------------------------------------------------
// This utility function implements a simple model of Torsen limited-slip
// differential with a max_bias:1 torque bias ratio.
// We hard-code the speed difference range over which the torque bias grows from
// a value of 1 to a value of max_bias to the interval [0.25, 0.5].
// -----------------------------------------------------------------------------
void differentialSplitXWD(double torque,
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
// -----------------------------------------------------------------------------
void ChSimpleDrivelineXWD::Synchronize(double time, const DriverInputs& driver_inputs, double driveshaft_torque) {
    if (!m_connected)
        return;

    double alpha = 1.0 / double(m_shaft_left.size());

    // Set driveshaft speed (output to transmission).
    // The axle shafts spin opposite to the wheels and are driven by the negative of the wheel torque, so that a positive
    // driveshaft speed and torque correspond to forward motion and traction (as in the shafts-based drivelines). The
    // conical gears reduce the driveshaft speed by their gear ratio and amplify the torque by the inverse ratio.
    m_driveshaft_speed = 0;
    for (int i = 0; i < m_shaft_left.size(); i++) {
        m_driveshaft_speed -= alpha * 0.5 * (m_shaft_left[i]->GetPosDt() + m_shaft_right[i]->GetPosDt());
    }
    m_driveshaft_speed /= GetConicalGearRatio();

    // Split the input torque over all driven axles: equally, or (if the inter-axle distribution is locked) with a bias
    // that grows with the axle speed excess over the slowest axle, as for a Torsen differential.
    size_t num_axles = m_shaft_left.size();
    std::vector<double> axle_fraction(num_axles, alpha);
    if (m_central_locked && num_axles > 1) {
        std::vector<double> speed(num_axles);
        double min_speed = 1e30;
        for (size_t i = 0; i < num_axles; i++) {
            speed[i] = std::abs(0.5 * (m_shaft_left[i]->GetPosDt() + m_shaft_right[i]->GetPosDt()));
            min_speed = std::min(min_speed, speed[i]);
        }
        double sum = 0;
        for (size_t i = 0; i < num_axles; i++) {
            double diff = speed[i] - min_speed;
            double bias = 1;
            if (diff > 0.5)
                bias = kLockedBiasXWD;
            else if (diff > 0.25)
                bias = 4 * (kLockedBiasXWD - 1) * diff + (2 - kLockedBiasXWD);
            axle_fraction[i] = 1 / bias;
            sum += axle_fraction[i];
        }
        for (auto& f : axle_fraction)
            f /= sum;
    }

    // Split the axle torques for the corresponding left/right wheels and apply
    // them to the suspension wheel shafts.

    for (int axle = 0; axle < m_shaft_left.size(); axle++) {
        double torque_left;
        double torque_right;
        double torque_axle = axle_fraction[axle] * driveshaft_torque / GetConicalGearRatio();
        bool locked = m_all_axles_locked || (axle < (int)m_axle_locked.size() && m_axle_locked[axle]);
        differentialSplitXWD(torque_axle, locked ? kLockedBiasXWD : GetDifferentialMaxBias(),
                             m_shaft_left[axle]->GetPosDt(), m_shaft_right[axle]->GetPosDt(), torque_left,
                             torque_right);
        m_shaft_left[axle]->SetAppliedLoad(-torque_left);
        m_shaft_right[axle]->SetAppliedLoad(-torque_right);
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
double ChSimpleDrivelineXWD::GetSpindleTorque(int axle, VehicleSide side) const {
    if (!m_connected)
        return 0;

    switch (side) {
        case LEFT:
            return -m_shaft_left[axle]->GetAppliedLoad();
        case RIGHT:
            return -m_shaft_right[axle]->GetAppliedLoad();
    }

    return 0;
}

// -----------------------------------------------------------------------------
void ChSimpleDrivelineXWD::Disconnect() {
    m_connected = false;
}

}  // end namespace vehicle
}  // end namespace chrono
