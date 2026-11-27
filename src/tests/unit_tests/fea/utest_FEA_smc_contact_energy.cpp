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
// Test for SMC contacts of FEA contact surfaces: no energy is created by the contact.
//
// A soft FE part (ChContactSurfaceMesh on its boundary) is dropped onto a fixed box with a
// restitution coefficient 0.2. The damping of the contact force is computed from the mass of
// the contactables, i.e. the element mass lumped at the nodes of each contact triangle. The
// mechanical energy (kinetic + gravitational) after the impact must not exceed the energy at
// release, and the part must come to rest on the box.
// Parts: blocks of 4-node and 10-node tetrahedra and 20-node hexahedra (step 5e-4 s), and a plate of
// Reissner shell elements (step 1e-4 s: the thin plate has light nodes, and the explicit contact force
// needs a smaller step). The 10-node, 20-node and shell elements have no own nodal mass lumping:
// their masses are lumped from the element mass matrices.
//
// =============================================================================

#include <array>
#include <cmath>
#include <map>

#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/solver/ChDirectSolverLS.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChElementTetraCorot_4.h"
#include "chrono/fea/ChElementTetraCorot_10.h"
#include "chrono/fea/ChElementHexaCorot_20.h"
#include "chrono/fea/ChElementShellReissner4.h"
#include "chrono/fea/ChMaterialShellReissner.h"
#include "chrono/fea/ChContactSurfaceMesh.h"

#include "gtest/gtest.h"

using namespace chrono;
using namespace chrono::fea;

enum class PartType { TETRA4, TETRA10, HEXA20, REISSNER };

static const double rho = 1000;  // density
static const double z0 = 0.25;   // release height of the bottom of the part

// Block 0.3 x 0.08 x 0.08 m of solid elements, built on a regular grid of hexahedral cells (6 tetrahedra per cell)
static double CreateBlock(std::shared_ptr<ChMesh> mesh, PartType type) {
    auto material = chrono_types::make_shared<ChContinuumElastic>();
    material->SetYoungModulus(2e6);
    material->SetPoissonRatio(0.3);
    material->SetDensity(rho);
    material->SetRayleighDampingBeta(2e-3);

    const int nx = 6, ny = 2, nz = 2;
    const double L[3] = {0.3, 0.08, 0.08};
    const double d[3] = {L[0] / nx, L[1] / ny, L[2] / nz};

    // Nodes on a grid of half cells (corners and mid-edge nodes), created on demand
    std::map<std::array<int, 3>, std::shared_ptr<ChNodeFEAxyz>> grid;
    auto node = [&](const std::array<int, 3>& g) {
        auto it = grid.find(g);
        if (it != grid.end())
            return it->second;
        auto n = chrono_types::make_shared<ChNodeFEAxyz>(ChVector3d(-0.15 + g[0] * d[0] / 2, -0.04 + g[1] * d[1] / 2, z0 + g[2] * d[2] / 2));
        grid[g] = n;
        mesh->AddNode(n);
        return n;
    };
    auto mid = [&](const std::array<int, 3>& a, const std::array<int, 3>& b) { return node({(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2}); };

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                // cell corners, in the order of the hexahedron elements
                std::array<int, 3> c[8] = {{2 * i, 2 * j, 2 * k},     {2 * i + 2, 2 * j, 2 * k},     {2 * i + 2, 2 * j + 2, 2 * k},     {2 * i, 2 * j + 2, 2 * k},
                                           {2 * i, 2 * j, 2 * k + 2}, {2 * i + 2, 2 * j, 2 * k + 2}, {2 * i + 2, 2 * j + 2, 2 * k + 2}, {2 * i, 2 * j + 2, 2 * k + 2}};
                if (type == PartType::HEXA20) {
                    const int edges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
                    std::shared_ptr<ChNodeFEAxyz> n[20];
                    for (int q = 0; q < 8; q++)
                        n[q] = node(c[q]);
                    for (int q = 0; q < 12; q++)
                        n[8 + q] = mid(c[edges[q][0]], c[edges[q][1]]);
                    auto element = chrono_types::make_shared<ChElementHexaCorot_20>();
                    element->SetNodes(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10], n[11], n[12], n[13], n[14], n[15], n[16], n[17], n[18], n[19]);
                    element->SetMaterial(material);
                    mesh->AddElement(element);
                    continue;
                }
                std::array<int, 3> cb[8] = {c[0], c[1], c[3], c[2], c[4], c[5], c[7], c[6]};
                const int tets[6][4] = {{0, 1, 3, 7}, {0, 1, 7, 5}, {0, 5, 7, 4}, {0, 3, 2, 7}, {0, 2, 6, 7}, {0, 6, 4, 7}};
                for (const auto& t : tets) {
                    std::array<int, 3> p[4] = {cb[t[0]], cb[t[1]], cb[t[2]], cb[t[3]]};
                    ChVector3d P0 = node(p[0])->GetPos(), P1 = node(p[1])->GetPos();
                    ChVector3d P2 = node(p[2])->GetPos(), P3 = node(p[3])->GetPos();
                    if (Vdot(Vcross(P1 - P0, P2 - P0), P3 - P0) < 0)
                        std::swap(p[1], p[2]);
                    if (type == PartType::TETRA4) {
                        auto element = chrono_types::make_shared<ChElementTetraCorot_4>();
                        element->SetNodes(node(p[0]), node(p[1]), node(p[2]), node(p[3]));
                        element->SetMaterial(material);
                        mesh->AddElement(element);
                    } else {
                        auto element = chrono_types::make_shared<ChElementTetraCorot_10>();
                        element->SetNodes(node(p[0]), node(p[1]), node(p[2]), node(p[3]), mid(p[0], p[1]), mid(p[1], p[2]), mid(p[2], p[0]), mid(p[0], p[3]), mid(p[1], p[3]),
                                          mid(p[2], p[3]));
                        element->SetMaterial(material);
                        mesh->AddElement(element);
                    }
                }
            }
        }
    }
    return rho * L[0] * L[1] * L[2];
}

// Square plate 0.3 x 0.3 m, thickness 5 mm, of Reissner shells
static double CreatePlate(std::shared_ptr<ChMesh> mesh) {
    const double L = 0.3, t = 0.005;
    {
        const int n = 6;
        auto material = chrono_types::make_shared<ChMaterialShellReissnerIsothropic>(rho, 2e6, 0.3);
        std::vector<std::shared_ptr<ChNodeFEAxyzrot>> nodes;
        for (int j = 0; j <= n; j++)
            for (int i = 0; i <= n; i++) {
                auto node = chrono_types::make_shared<ChNodeFEAxyzrot>(ChFrame<>(ChVector3d(-L / 2 + i * L / n, -L / 2 + j * L / n, z0)));
                mesh->AddNode(node);
                nodes.push_back(node);
            }
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                auto element = chrono_types::make_shared<ChElementShellReissner4>();
                mesh->AddElement(element);
                element->SetNodes(nodes[j * (n + 1) + i], nodes[j * (n + 1) + i + 1], nodes[(j + 1) * (n + 1) + i + 1], nodes[(j + 1) * (n + 1) + i]);
                element->AddLayer(t, 0, material);
            }
    }
    return rho * L * L * t;
}

class FEA_SMC_contact : public ::testing::TestWithParam<PartType> {};

TEST_P(FEA_SMC_contact, no_energy_gain) {
    const double g = 9.81;

    ChSystemSMC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, -g));
    sys.SetContactForceModel(ChSystemSMC::Hertz);
    sys.SetSolver(chrono_types::make_shared<ChSolverSparseLU>());
    sys.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);

    auto ground_mat = chrono_types::make_shared<ChContactMaterialSMC>();
    ground_mat->SetYoungModulus(1e7f);
    ground_mat->SetFriction(0.3f);
    ground_mat->SetRestitution(0.5f);
    auto ground = chrono_types::make_shared<ChBodyEasyBox>(2, 2, 0.2, 1000, false, true, ground_mat);
    ground->SetPos(ChVector3d(0, 0, -0.1));
    ground->SetFixed(true);
    sys.Add(ground);

    auto mesh = chrono_types::make_shared<ChMesh>();
    PartType type = GetParam();
    bool block = (type == PartType::TETRA4 || type == PartType::TETRA10 || type == PartType::HEXA20);
    double total_mass = block ? CreateBlock(mesh, type) : CreatePlate(mesh);
    double h = block ? 5e-4 : 1e-4;

    auto contact_mat = chrono_types::make_shared<ChContactMaterialSMC>();
    contact_mat->SetYoungModulus(1e7f);
    contact_mat->SetFriction(0.3f);
    contact_mat->SetRestitution(0.2f);
    auto surface = chrono_types::make_shared<ChContactSurfaceMesh>(contact_mat);
    surface->AddFacesFromBoundary(*mesh, 0.002);
    mesh->AddContactSurface(surface);
    sys.Add(mesh);

    // Mechanical energy, with the mass distributed equally over the nodes
    const auto& nodes = mesh->GetNodes();
    auto position = [](const std::shared_ptr<ChNodeFEAbase>& node, ChVector3d& vel) {
        if (auto n = std::dynamic_pointer_cast<ChNodeFEAxyz>(node)) {
            vel = n->GetPosDt();
            return n->GetPos();
        }
        auto n = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(node);
        vel = n->GetPosDt();
        return n->GetPos();
    };
    auto energy = [&](double& z_mean) {
        double E = 0;
        double m = total_mass / nodes.size();
        z_mean = 0;
        for (const auto& node : nodes) {
            ChVector3d vel;
            ChVector3d pos = position(node, vel);
            E += 0.5 * m * vel.Length2() + m * g * pos.z();
            z_mean += pos.z() / nodes.size();
        }
        return E;
    };
    double z_mean;
    double E0 = energy(z_mean);

    bool touched = false;
    double E_max = -1e10;
    while (sys.GetChTime() < 0.8) {
        sys.DoStepDynamics(h);
        if (sys.GetNumContacts() > 0)
            touched = true;
        if (touched)
            E_max = std::max(E_max, energy(z_mean));
    }
    energy(z_mean);

    // Lumped nodal masses: positive, summing to the element masses (HexaCorot_20 lumps its own mass matrix,
    // which uses its numerically integrated volume)
    double mass_sum = 0;
    for (const auto& node : nodes) {
        EXPECT_GT(node->m_TotalMass, 0);
        mass_sum += node->m_TotalMass;
    }
    if (type != PartType::HEXA20)
        EXPECT_NEAR(mass_sum, total_mass, 1e-9 * total_mass);

    std::cout << "energy at release " << E0 << "  max energy after impact " << E_max << "  final mean z " << z_mean << "  lumped mass " << mass_sum << " (" << total_mass << ")"
              << std::endl;
    ASSERT_TRUE(touched);
    ASSERT_LT(E_max, E0);
    double z_rest = block ? 0.04 : 0.0;
    ASSERT_NEAR(z_mean, z_rest, 0.006);
}

INSTANTIATE_TEST_SUITE_P(FEA, FEA_SMC_contact, ::testing::Values(PartType::TETRA4, PartType::TETRA10, PartType::HEXA20, PartType::REISSNER));
