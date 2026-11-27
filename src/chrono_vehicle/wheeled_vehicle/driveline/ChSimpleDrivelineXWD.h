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

#ifndef CH_SIMPLE_DRIVELINE_XWD_H
#define CH_SIMPLE_DRIVELINE_XWD_H

#include "chrono_vehicle/ChApiVehicle.h"
#include "chrono_vehicle/wheeled_vehicle/ChDrivelineWV.h"

namespace chrono {
namespace vehicle {

/// @addtogroup vehicle_wheeled_driveline
/// @{

/// Simple driveline model. This template can be used to model a XWD driveline.
/// Number of axles can be 1 to X.
/// It uses a constant torque split depending on the number of axles driven and a
/// simple model for Torsen limited-slip differentials.
class CH_VEHICLE_API ChSimpleDrivelineXWD : public ChDrivelineWV {
  public:
    ChSimpleDrivelineXWD(const std::string& name);

    virtual ~ChSimpleDrivelineXWD() {}

    /// Get the name of the vehicle subsystem template.
    virtual std::string GetTemplateName() const override { return "SimpleDrivelineXWD"; }

    /// Lock/unlock the differential on the specified driven axle (-1 for all axles).
    /// A locked differential is modeled as a limited-slip differential with a very large torque bias ratio, so that
    /// the torque goes to the slower wheel as soon as the two wheel speeds differ.
    virtual void LockAxleDifferential(int axle, bool lock) override;

    /// Lock/unlock the inter-axle torque distribution (the 'which' argument is ignored).
    /// When locked, the driveshaft torque is distributed with a very large bias towards the slowest axles instead of
    /// equally over all driven axles.
    virtual void LockCentralDifferential(int which, bool lock) override;

    /// Return the number of driven axles.
    virtual unsigned int GetNumDrivenAxles() const override { return static_cast<int>(m_shaft_left.size()); }

    /// Initialize the driveline subsystem.
    /// This function connects this driveline subsystem to the specified axle subsystems.
    virtual void Initialize(std::shared_ptr<ChChassis> chassis,   ///< associated chassis subsystem
                            const ChAxleList& axles,              ///< list of all vehicle axle subsystems
                            const std::vector<int>& driven_axles  ///< indexes of the driven vehicle axles
                            ) override;

    /// Update the driveline subsystem: apply the specified motor torque.
    /// This represents the input to the driveline subsystem from the powertrain system.
    virtual void Synchronize(double time,                        ///< current time
                             const DriverInputs& driver_inputs,  ///< current driver inputs
                             double driveshaft_torque            ///< input transmission torque
                             ) override;

    /// Get the motor torque to be applied to the specified spindle.
    virtual double GetSpindleTorque(int axle, VehicleSide side) const override;

    /// Disconnect driveline from driven wheels.
    virtual void Disconnect() override;

    /// Return the output driveline speed of the driveshaft.
    /// This represents the output from the driveline subsystem that is passed to the transmission subsystem.
    virtual double GetOutputDriveshaftSpeed() const override { return m_driveshaft_speed; }

  protected:
    /// Return the torque bias ratio every axle-wise differential.
    /// This is a simple model of a Torsen limited-slip differential.
    virtual double GetDifferentialMaxBias() const = 0;

    /// Return the gear ratio of the conical gear of each driven axle (differential speed / driveshaft speed).
    /// The default value of 1 corresponds to no reduction between driveshaft and axles.
    virtual double GetConicalGearRatio() const { return 1; }

  private:
    bool m_connected;
    std::vector<bool> m_axle_locked;                       ///< differential locked, per driven axle
    bool m_all_axles_locked;                               ///< all axle differentials locked
    bool m_central_locked;                                 ///< inter-axle distribution locked
    double m_driveshaft_speed;                             ///< output to transmission
    std::vector<std::shared_ptr<ChShaft> > m_shaft_left;   ///< associated left axles
    std::vector<std::shared_ptr<ChShaft> > m_shaft_right;  ///< associated right axles
};

/// @} vehicle_wheeled_driveline

}  // end namespace vehicle
}  // end namespace chrono

#endif
