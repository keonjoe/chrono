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
//
// Unit test for transient heat conduction with the Newton-iterated implicit
// integrators, on a bar of hexahedra with the temperature imposed at x = 0.
//
// - Prediction: with the temperature also imposed at x = L, the linear steady
//   profile must be kept by the prediction of the new state alone (no Newton
//   iteration). The temperature is stored at the velocity level of the field,
//   so the prediction must include the conduction flux through dF/dv.
// - Convection: with a convecting end at x = L, the steady end temperature is
//      T_L = T_inf + (T_0 - T_inf) / (1 + h*L/k)
//   (exact for linear elements), independent of the step size, and the Newton
//   iteration converges in few iterations even for h*A*dt much larger than the
//   heat capacity of the end nodes, as the convection load provides its Jacobian.
// - Radiation: with a radiating end at x = L exchanging heat with an environment
//   at T_env, the steady end temperature solves
//      k*(T_L - T_0)/L = e*sigma*(T_env^4 - T_L^4)
//   here with T_env > T_0, so that the end is heated by the environment.
// - Transient: Newmark and HHT start from the rate of change of temperature
//   stored in the field (zero at the start), so two identical systems give the
//   same response, which agrees with the series solution of the convecting bar.
//
// =============================================================================

#include <cmath>
#include <string>

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/solver/ChDirectSolverLS.h"
#include "chrono/timestepper/ChTimestepperImplicit.h"
#include "chrono/timestepper/ChTimestepperHHT.h"
#include "chrono/fea/multiphysics/ChFEModelThermal.h"
#include "chrono/fea/multiphysics/ChMaterial3DThermalLinear.h"
#include "chrono/fea/multiphysics/ChBuilderVolume.h"
#include "chrono/fea/multiphysics/ChFieldElementLoadableSurface.h"
#include "chrono/fea/multiphysics/ChLoaderHeatConvection.h"
#include "chrono/fea/multiphysics/ChLoaderHeatRadiation.h"

#include "gtest/gtest.h"

using namespace chrono;
using namespace chrono::fea;

// Bar 0.2 m long, 0.02 x 0.02 m section, k = 50 W/(m K), rho*c = 1e5 J/(m^3 K), convection h = 500 W/(m^2 K).
// With 20 elements, the heat capacity of the end nodes is h*A*1s: dt = 50 s gives h*A*dt/C = 50.
const double L = 0.2;
const double W = 0.02;
const double k_cond = 50;
const double rho = 1000;
const double c_heat = 100;
const double h_conv = 500;
const double T_0 = 400;
const double T_inf = 300;
const int nx = 20;
const double sigma_SB = 5.670374419e-8;

enum class EndCondition { FIXED, CONVECTION, RADIATION };

// End temperature of the bar with a convecting end at time t, for the uniform initial temperature T_0: series of the
// modes sin(mu*x/L), tan(mu) = -mu/Bi, Bi = h*L/k
double ConvectionSeries(double t) {
    const double Bi = h_conv * L / k_cond;
    const double alpha = k_cond / (rho * c_heat);
    const double T_L = T_inf + (T_0 - T_inf) / (1 + Bi);
    double sum = 0;
    for (int n = 1; n <= 200; n++) {
        double lo = (n - 0.5) * CH_PI, hi = n * CH_PI;
        for (int it = 0; it < 100; it++) {
            double mu = 0.5 * (lo + hi);
            if ((std::sin(lo) + lo / Bi * std::cos(lo)) * (std::sin(mu) + mu / Bi * std::cos(mu)) <= 0)
                hi = mu;
            else
                lo = mu;
        }
        double mu = 0.5 * (lo + hi);
        double b = ((std::sin(mu) - mu * std::cos(mu)) / (mu * mu)) / (0.5 - std::sin(2 * mu) / (4 * mu));
        sum += b * std::sin(mu) * std::exp(-alpha * mu * mu * t / (L * L));
    }
    return T_L + (T_0 - T_L) * sum;
}

class ThermalBar {
  public:
    ThermalBar(const std::string& stepper, EndCondition end, double T_root = T_0, double T_ext = T_inf) {
        bool fixed_end = (end == EndCondition::FIXED);
        field = chrono_types::make_shared<ChFieldTemperature>();
        sys.Add(field);
        auto model = chrono_types::make_shared<ChFEModelThermal>(field);
        sys.Add(model);
        auto material = chrono_types::make_shared<ChMaterial3DThermalLinear>();
        material->SetDensity(rho);
        material->SetSpecificHeatCapacity(c_heat);
        material->SetThermalConductivity(k_cond);
        model->material = material;

        builder.BuildVolume(ChFrame<>(), nx, 1, 1, L, W, W);
        builder.AddToModel(model);
        for (int i = 0; i <= nx; i++) {
            for (int j = 0; j <= 1; j++) {
                for (int l = 0; l <= 1; l++) {
                    auto& data = field->NodeData(builder.nodes.at(i, j, l));
                    if (fixed_end) {
                        // linear steady profile from T_root to T_ext
                        data.T() = T_root + (T_ext - T_root) * i / nx;
                        data.SetFixed(i == 0 || i == nx);
                    } else {
                        data.T() = T_root;
                        data.SetFixed(i == 0);
                    }
                }
            }
        }

        if (!fixed_end) {
            auto loads = chrono_types::make_shared<ChLoadContainer>();
            sys.Add(loads);
            for (auto& face : builder.faces_x_hi) {
                auto loadable = chrono_types::make_shared<ChFieldElementLoadableSurface>(face, field);
                if (end == EndCondition::CONVECTION) {
                    auto conv = chrono_types::make_shared<ChLoaderHeatConvection>(loadable, field);
                    conv->SetSurfaceConvectionCoeff(h_conv);
                    conv->SetFluidTemperature(T_ext);
                    conv->SetIntegrationPoints(2);
                    loads->Add(conv);
                } else {
                    auto rad = chrono_types::make_shared<ChLoaderHeatRadiation>(loadable, field);
                    rad->SetSurfaceEmissivity(1.0);
                    rad->SetEnvironmentTemperature(T_ext);
                    rad->SetIntegrationPoints(2);
                    loads->Add(rad);
                }
            }
        }

        sys.SetSolver(chrono_types::make_shared<ChSolverSparseQR>());
        if (stepper == "EULER_IMPLICIT") {
            auto t = chrono_types::make_shared<ChTimestepperEulerImplicit>(&sys);
            sys.SetTimestepper(t);
            integrator = t;
        } else if (stepper == "TRAPEZOIDAL") {
            auto t = chrono_types::make_shared<ChTimestepperTrapezoidal>(&sys);
            sys.SetTimestepper(t);
            integrator = t;
        } else if (stepper == "NEWMARK") {
            auto t = chrono_types::make_shared<ChTimestepperNewmark>(&sys);
            sys.SetTimestepper(t);
            integrator = t;
        } else {
            auto t = chrono_types::make_shared<ChTimestepperHHT>(&sys);
            t->SetStepControl(false);
            sys.SetTimestepper(t);
            integrator = t;
        }
    }

    double T(int i) { return field->NodeData(builder.nodes.at(i, 0, 0)).T(); }

    ChSystemNSC sys;
    std::shared_ptr<ChFieldTemperature> field;
    ChBuilderVolumeBox builder;
    std::shared_ptr<ChTimestepperImplicit> integrator;
};

// Steady linear profile between two imposed temperatures, kept by the prediction alone
class ThermalPredictionTest : public ::testing::TestWithParam<std::string> {};

TEST_P(ThermalPredictionTest, steady_profile) {
    ThermalBar bar(GetParam(), EndCondition::FIXED);
    bar.integrator->SetMaxIters(0);
    for (int step = 0; step < 5; step++)
        bar.sys.DoStepDynamics(10.0);
    for (int i = 0; i <= nx; i++)
        ASSERT_NEAR(bar.T(i), T_0 + (T_inf - T_0) * i / nx, 1e-9) << "node " << i;
}

INSTANTIATE_TEST_SUITE_P(FEA, ThermalPredictionTest, ::testing::Values("EULER_IMPLICIT", "TRAPEZOIDAL"));

// Steady end temperature of a bar with a convecting end, for large steps
class ThermalConvectionTest : public ::testing::TestWithParam<std::string> {};

TEST_P(ThermalConvectionTest, steady_end_temperature) {
    const double T_L = T_inf + (T_0 - T_inf) / (1 + h_conv * L / k_cond);  // 333.333 K
    const double dt = 50;
    const int num_steps = 30;  // 100 times the decay time of the slowest mode

    ThermalBar bar(GetParam(), EndCondition::CONVECTION);
    for (int step = 0; step < num_steps; step++)
        bar.sys.DoStepDynamics(dt);

    // (Newmark with gamma = 0.6 damps the stiff conduction modes slowly at this step)
    ASSERT_NEAR(bar.T(nx), T_L, GetParam() == "NEWMARK" ? 1e-5 : 1e-10);
    ASSERT_EQ(bar.integrator->GetNumTerminated(), 0u);
    ASSERT_LE(bar.integrator->GetNumIterations(), 3u * num_steps);
}

// (The trapezoidal rule is not L-stable: at such steps the stiff conduction modes decay too slowly for this check.)
INSTANTIATE_TEST_SUITE_P(FEA, ThermalConvectionTest, ::testing::Values("EULER_IMPLICIT", "NEWMARK", "HHT"));

// Steady end temperature of a bar heated by radiation from a hotter environment
class ThermalRadiationTest : public ::testing::TestWithParam<std::string> {};

TEST_P(ThermalRadiationTest, steady_end_temperature) {
    const double T_root = 300;
    const double T_env = 1000;
    double lo = T_root, hi = T_env;  // k*(T_L - T_root)/L - sigma*(T_env^4 - T_L^4) is increasing in T_L
    for (int it = 0; it < 200; it++) {
        double T = 0.5 * (lo + hi);
        if (k_cond * (T - T_root) / L - sigma_SB * (std::pow(T_env, 4) - std::pow(T, 4)) > 0)
            hi = T;
        else
            lo = T;
    }
    const double T_L = 0.5 * (lo + hi);  // 511.312 K

    ThermalBar bar(GetParam(), EndCondition::RADIATION, T_root, T_env);
    for (int step = 0; step < 30; step++)
        bar.sys.DoStepDynamics(50);

    ASSERT_NEAR(bar.T(nx), T_L, 1e-8);
    ASSERT_EQ(bar.integrator->GetNumTerminated(), 0u);
}

INSTANTIATE_TEST_SUITE_P(FEA, ThermalRadiationTest, ::testing::Values("EULER_IMPLICIT", "HHT"));

// The field nodes store the rate of change of temperature as acceleration-level state: zero at the start, then the
// value scattered by the integrator
TEST(ThermalAccelerationTest, gather_scatter) {
    ThermalBar bar("NEWMARK", EndCondition::CONVECTION);
    bar.sys.Setup();
    ChStateDelta a(bar.sys.GetNumCoordsVelLevel(), &bar.sys);
    a.setConstant(1e9);
    bar.sys.StateGatherAcceleration(a);
    ASSERT_EQ(a.lpNorm<Eigen::Infinity>(), 0.0);

    for (int step = 0; step < 3; step++)
        bar.sys.DoStepDynamics(0.05);
    a.setConstant(1e9);
    bar.sys.StateGatherAcceleration(a);
    const auto& a_integrator = std::dynamic_pointer_cast<ChTimestepperIIorder>(bar.integrator)->GetStateAcc();
    ASSERT_EQ((a - a_integrator).lpNorm<Eigen::Infinity>(), 0.0);
    ASSERT_GT(a.lpNorm<Eigen::Infinity>(), 0.0);
}

// Transient of the convecting bar with Newmark and HHT: reproducible, and close to the series solution
class ThermalTransientTest : public ::testing::TestWithParam<std::string> {};

TEST_P(ThermalTransientTest, end_temperature) {
    const double dt = 0.05;
    const double t_end = 5;
    ThermalBar bar1(GetParam(), EndCondition::CONVECTION);
    ThermalBar bar2(GetParam(), EndCondition::CONVECTION);
    for (int step = 0; step < (int)std::lround(t_end / dt); step++) {
        bar1.sys.DoStepDynamics(dt);
        bar2.sys.DoStepDynamics(dt);
    }
    ASSERT_NEAR(bar1.T(nx), bar2.T(nx), 1e-9);
    ASSERT_NEAR(bar1.T(nx), ConvectionSeries(t_end), 0.2);
}

INSTANTIATE_TEST_SUITE_P(FEA, ThermalTransientTest, ::testing::Values("NEWMARK", "HHT"));
