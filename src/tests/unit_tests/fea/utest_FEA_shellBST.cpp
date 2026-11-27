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
// Tests for the BST thin shell triangle (ChElementShellBST).
//
// - Curvature patch test: a quadratic deflection w(x,y) imposed on all nodes of a flat mesh must give its constant
//   curvature in every element that has three neighbours, on meshes with parallel and with alternating diagonals.
// - Linear static plate benchmarks against Kirchhoff plate theory (Timoshenko & Woinowsky-Krieger):
//     simply supported square plate, uniform pressure q:  w_max = 0.00406235 q a^4 / D
//     clamped square plate, uniform pressure q:           w_max = 0.00126532 q a^4 / D
//     clamped circular plate of radius R:                 w_max = q R^4 / (64 D)
//     cantilever strip of length L, width b, tip load P:  w_max = P L^3 / (3 D b)   (nu = 0)
//   with D = E t^3 / (12 (1 - nu^2)). The sum of the support reactions must balance the applied load.
//
// =============================================================================

#include <array>
#include <cmath>
#include <map>
#include <vector>

#include "chrono/fea/ChElementShellBST.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChNodeFEAxyz.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/solver/ChDirectSolverLS.h"

#include "gtest/gtest.h"

using namespace chrono;
using namespace chrono::fea;

namespace {

const double E = 2.1e11;
const double t = 0.002;

struct TriMesh {
    std::vector<ChVector3d> p;
    std::vector<std::array<int, 3>> tri;
};

// Rectangle [0,a]x[0,b] with nx x ny cells; each cell split along the diagonal (0,0)-(1,1), or along alternating
// diagonals ("union jack") if alternating = true.
TriMesh RectangleMesh(double a, double b, int nx, int ny, bool alternating) {
    TriMesh m;
    for (int j = 0; j <= ny; j++)
        for (int i = 0; i <= nx; i++)
            m.p.push_back(ChVector3d(a * i / nx, b * j / ny, 0));
    auto id = [&](int i, int j) { return j * (nx + 1) + i; };
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            int n0 = id(i, j), n1 = id(i + 1, j), n2 = id(i + 1, j + 1), n3 = id(i, j + 1);
            if (alternating && (i + j) % 2 == 1) {
                m.tri.push_back({n0, n1, n3});
                m.tri.push_back({n1, n2, n3});
            } else {
                m.tri.push_back({n0, n1, n2});
                m.tri.push_back({n0, n2, n3});
            }
        }
    return m;
}

// Disk of radius R with nr rings of nodes (6 k nodes on ring k).
TriMesh DiskMesh(double R, int nr) {
    TriMesh m;
    std::vector<std::vector<int>> ring(nr + 1);
    m.p.push_back(VNULL);
    ring[0] = {0};
    for (int k = 1; k <= nr; k++)
        for (int s = 0; s < 6 * k; s++) {
            double a = CH_2PI * s / (6 * k);
            ring[k].push_back((int)m.p.size());
            m.p.push_back(ChVector3d(R * k / nr * std::cos(a), R * k / nr * std::sin(a), 0));
        }
    for (int s = 0; s < 6; s++)
        m.tri.push_back({0, ring[1][s], ring[1][(s + 1) % 6]});
    for (int k = 2; k <= nr; k++) {
        int ni = 6 * (k - 1), no = 6 * k;
        int a = 0, b = 0;
        while (a < ni || b < no) {
            if (b < no && (a >= ni || (b + 1.0) / no <= (a + 1.0) / ni + 1e-12)) {
                m.tri.push_back({ring[k - 1][a % ni], ring[k][b], ring[k][(b + 1) % no]});
                b++;
            } else {
                m.tri.push_back({ring[k - 1][a % ni], ring[k][b % no], ring[k - 1][(a + 1) % ni]});
                a++;
            }
        }
    }
    return m;
}

// FEA model of a triangle mesh, with the neighbouring nodes of each element found from the shared edges.
struct ShellModel {
    ChSystemNSC sys;
    std::shared_ptr<ChMesh> mesh;
    std::vector<std::shared_ptr<ChNodeFEAxyz>> nodes;
    std::vector<std::shared_ptr<ChElementShellBST>> elements;
    std::vector<bool> on_boundary;

    // clamped(p, q): true if the boundary edge p-q is clamped
    template <typename Clamped>
    ShellModel(const TriMesh& m, double nu, Clamped clamped) {
        mesh = chrono_types::make_shared<ChMesh>();
        sys.Add(mesh);
        auto elasticity = chrono_types::make_shared<ChElasticityKirchhoffIsothropic>(E, nu);
        auto material = chrono_types::make_shared<ChMaterialShellKirchhoff>(elasticity);
        material->SetDensity(7800);
        for (auto& p : m.p) {
            nodes.push_back(chrono_types::make_shared<ChNodeFEAxyz>(p));
            mesh->AddNode(nodes.back());
        }
        // edge (i,j) -> (element, opposite node)
        std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> edges;
        for (int e = 0; e < (int)m.tri.size(); e++)
            for (int k = 0; k < 3; k++) {
                int i = m.tri[e][(k + 1) % 3], j = m.tri[e][(k + 2) % 3];
                edges[{std::min(i, j), std::max(i, j)}].push_back({e, m.tri[e][k]});
            }
        on_boundary.assign(m.p.size(), false);
        for (auto& ed : edges)
            if (ed.second.size() == 1)
                on_boundary[ed.first.first] = on_boundary[ed.first.second] = true;
        for (int e = 0; e < (int)m.tri.size(); e++) {
            std::shared_ptr<ChNodeFEAxyz> nb[3];
            for (int k = 0; k < 3; k++) {
                int i = m.tri[e][(k + 1) % 3], j = m.tri[e][(k + 2) % 3];
                for (auto& o : edges[{std::min(i, j), std::max(i, j)}])
                    if (o.first != e)
                        nb[k] = nodes[o.second];
            }
            auto element = chrono_types::make_shared<ChElementShellBST>();
            element->SetNodes(nodes[m.tri[e][0]], nodes[m.tri[e][1]], nodes[m.tri[e][2]], nb[0], nb[1], nb[2]);
            element->AddLayer(t, 0, material);
            for (int k = 0; k < 3; k++)
                if (!nb[k] && clamped(m.p[m.tri[e][(k + 1) % 3]], m.p[m.tri[e][(k + 2) % 3]]))
                    element->SetEdgeClamped(k, true);
            mesh->AddElement(element);
            elements.push_back(element);
        }
        sys.SetSolver(chrono_types::make_shared<ChSolverSparseLU>());
    }

    // Linear static solution for the nodal loads F; returns the sum of the reactions at the fixed nodes.
    ChVector3d SolveStatic(const std::vector<ChVector3d>& F) {
        for (size_t i = 0; i < nodes.size(); i++)
            nodes[i]->SetForce(F[i]);
        sys.DoStaticLinear();
        // reaction at a fixed node: minus the element forces and the external load on it
        std::map<ChNodeFEAbase*, ChVector3d> f_int;
        for (auto& element : elements) {
            ChVectorDynamic<> Fi(element->GetNumCoordsPosLevel());
            element->ComputeInternalForces(Fi);
            for (unsigned int k = 0; k < element->GetNumNodes(); k++)
                f_int[element->GetNode(k).get()] += ChVector3d(Fi.segment(3 * k, 3));
        }
        ChVector3d R = VNULL;
        for (size_t i = 0; i < nodes.size(); i++)
            if (nodes[i]->IsFixed())
                R -= f_int[nodes[i].get()] + F[i];
        return R;
    }

    // Lumped nodal loads of a uniform pressure q along -z.
    std::vector<ChVector3d> PressureLoads(const TriMesh& m, double q) const {
        std::vector<ChVector3d> F(m.p.size(), VNULL);
        for (auto& tr : m.tri) {
            double A = 0.5 * Vcross(m.p[tr[1]] - m.p[tr[0]], m.p[tr[2]] - m.p[tr[0]]).Length();
            for (int k = 0; k < 3; k++)
                F[tr[k]].z() -= q * A / 3;
        }
        return F;
    }

    double MaxDeflection() const {
        double w = 0;
        for (auto& n : nodes)
            w = std::max(w, std::abs(n->GetPos().z() - n->GetX0().z()));
        return w;
    }
};

auto free_edges = [](const ChVector3d&, const ChVector3d&) { return false; };
auto clamped_edges = [](const ChVector3d&, const ChVector3d&) { return true; };

double PlateStiffness(double nu) {
    return E * t * t * t / (12 * (1 - nu * nu));
}

}  // end namespace

// -----------------------------------------------------------------------------

TEST(ShellBST, curvature_patch_test) {
    for (bool alternating : {false, true}) {
        TriMesh m = RectangleMesh(1.0, 1.0, 6, 6, alternating);
        ShellModel model(m, 0.3, free_edges);
        // initialize the elements in the flat reference configuration (all nodes fixed but one)
        for (size_t i = 1; i < model.nodes.size(); i++)
            model.nodes[i]->SetFixed(true);
        model.sys.DoStaticLinear();

        // curvature (w_xx, w_yy, 2 w_xy) of the imposed fields, in global axes
        const ChVector3d fields[3] = {ChVector3d(1, 0, 0), ChVector3d(0, 1, 0), ChVector3d(0, 0, 2)};
        const double eps = 1e-6;  // small deflection: linear kinematics
        for (int f = 0; f < 3; f++) {
            for (size_t i = 0; i < model.nodes.size(); i++) {
                double x = m.p[i].x(), y = m.p[i].y();
                double w = f == 0 ? 0.5 * x * x : (f == 1 ? 0.5 * y * y : x * y);
                model.nodes[i]->SetPos(m.p[i] + ChVector3d(0, 0, eps * w));
            }
            int checked = 0;
            for (auto& element : model.elements) {
                if (element->GetNumNodes() < 6)
                    continue;
                ChVectorDynamic<> Fi(element->GetNumCoordsPosLevel());
                element->ComputeInternalForces(Fi);
                // element curvature in global axes (element x axis along the edge from node 0 to node 1)
                ChVector3d x1 = (element->GetNodeMainTriangle(1)->GetX0() - element->GetNodeMainTriangle(0)->GetX0()).GetNormalized();
                double c = x1.x(), s = x1.y();
                double k11 = element->k.x() / eps, k22 = element->k.y() / eps, k12 = 0.5 * element->k.z() / eps;
                double g11 = c * c * k11 - 2 * c * s * k12 + s * s * k22;
                double g22 = s * s * k11 + 2 * c * s * k12 + c * c * k22;
                double g12 = c * s * k11 + (c * c - s * s) * k12 - c * s * k22;
                ASSERT_NEAR(g11, fields[f].x(), 1e-4) << "field " << f << " alternating " << alternating;
                ASSERT_NEAR(g22, fields[f].y(), 1e-4) << "field " << f << " alternating " << alternating;
                ASSERT_NEAR(2 * g12, fields[f].z(), 1e-4) << "field " << f << " alternating " << alternating;
                checked++;
            }
            ASSERT_GT(checked, 0);
        }
    }
}

TEST(ShellBST, boundary_fan) {
    // fan of triangles from one inner node to five nodes on the straight boundary y = 0: the triangles in the middle
    // of the fan have both neighbours' opposite vertices on the boundary, so their two edge bendings measure the
    // same curvature along the boundary
    double nu = 0.3;
    TriMesh m;
    for (int i = 0; i <= 4; i++)
        m.p.push_back(ChVector3d(0.25 * i, 0, 0));  // 0..4: boundary
    m.p.push_back(ChVector3d(0.5, 0.4, 0));         // 5: inner node
    m.p.push_back(ChVector3d(0, 1, 0));             // 6
    m.p.push_back(ChVector3d(1, 1, 0));             // 7
    for (int i = 0; i < 4; i++)
        m.tri.push_back({i, i + 1, 5});
    m.tri.push_back({0, 5, 6});
    m.tri.push_back({5, 4, 7});
    m.tri.push_back({5, 7, 6});
    ShellModel model(m, nu, free_edges);
    for (size_t i = 1; i < model.nodes.size(); i++)
        model.nodes[i]->SetFixed(true);
    model.sys.DoStaticLinear();

    const double eps = 1e-6;
    for (size_t i = 0; i < model.nodes.size(); i++)
        model.nodes[i]->SetPos(m.p[i] + ChVector3d(0, 0, eps * 0.5 * m.p[i].x() * m.p[i].x()));
    for (int e : {1, 2}) {
        auto& element = model.elements[e];
        ChVectorDynamic<> Fi(element->GetNumCoordsPosLevel());
        element->ComputeInternalForces(Fi);
        // element x axis along the edge from node 0 to node 1, i.e. along the boundary
        ASSERT_NEAR(element->k.x() / eps, 1.0, 1e-4);
        // minimum bending energy: no bending moment normal to the free boundary
        ASSERT_NEAR(element->k.y() / eps, -nu, 1e-4);
        ASSERT_NEAR(element->k.z() / eps, 0.0, 1e-4);
        ASSERT_TRUE(Fi.allFinite());
    }
}

TEST(ShellBST, simply_supported_square_plate) {
    double a = 0.5, q = 1e-2, nu = 0.3;
    double exact = 0.00406235 * q * std::pow(a, 4) / PlateStiffness(nu);
    double err_prev = 1;
    for (int n : {8, 16}) {
        TriMesh m = RectangleMesh(a, a, n, n, false);
        ShellModel model(m, nu, free_edges);
        for (size_t i = 0; i < model.nodes.size(); i++)
            model.nodes[i]->SetFixed(model.on_boundary[i]);
        auto F = model.PressureLoads(m, q);
        ChVector3d R = model.SolveStatic(F);
        double err = std::abs(model.MaxDeflection() / exact - 1);
        std::cout << "simply supported plate, " << n << "x" << n << " cells: w/w_exact = " << model.MaxDeflection() / exact << ", reaction / load = " << R.z() / (q * a * a)
                  << std::endl;
        ASSERT_LT(err, n == 8 ? 0.04 : 0.01);
        ASSERT_LT(err, err_prev);
        err_prev = err;
        ASSERT_NEAR(R.z(), q * a * a, 1e-6 * q * a * a);
    }
}

TEST(ShellBST, clamped_square_plate) {
    double a = 0.5, q = 1e-2, nu = 0.3;
    double exact = 0.00126532 * q * std::pow(a, 4) / PlateStiffness(nu);
    TriMesh m = RectangleMesh(a, a, 16, 16, false);
    ShellModel model(m, nu, clamped_edges);
    for (size_t i = 0; i < model.nodes.size(); i++)
        model.nodes[i]->SetFixed(model.on_boundary[i]);
    ChVector3d R = model.SolveStatic(model.PressureLoads(m, q));
    std::cout << "clamped plate: w/w_exact = " << model.MaxDeflection() / exact << std::endl;
    ASSERT_NEAR(model.MaxDeflection() / exact, 1.0, 0.06);
    ASSERT_NEAR(R.z(), q * a * a, 1e-6 * q * a * a);
}

TEST(ShellBST, clamped_circular_plate) {
    double radius = 0.25, q = 1e-2, nu = 0.3;
    double exact = q * std::pow(radius, 4) / (64 * PlateStiffness(nu));
    TriMesh m = DiskMesh(radius, 16);
    ShellModel model(m, nu, clamped_edges);
    for (size_t i = 0; i < model.nodes.size(); i++)
        model.nodes[i]->SetFixed(model.on_boundary[i]);
    auto F = model.PressureLoads(m, q);
    double load = 0;
    for (auto& f : F)
        load -= f.z();
    ChVector3d R = model.SolveStatic(F);
    std::cout << "clamped circular plate: w/w_exact = " << model.MaxDeflection() / exact << std::endl;
    ASSERT_NEAR(model.MaxDeflection() / exact, 1.0, 0.04);
    ASSERT_NEAR(R.z(), load, 1e-6 * load);
}

TEST(ShellBST, cantilever_strip) {
    double L = 0.5, b = 0.1, P = 1e-4;
    double exact = P * L * L * L / (3 * PlateStiffness(0) * b);
    int nx = 32, ny = 6;
    TriMesh m = RectangleMesh(L, b, nx, ny, false);
    auto clamped_root = [](const ChVector3d& p, const ChVector3d& q) { return p.x() < 1e-9 && q.x() < 1e-9; };
    ShellModel model(m, 0.0, clamped_root);
    std::vector<ChVector3d> F(m.p.size(), VNULL);
    for (size_t i = 0; i < m.p.size(); i++) {
        model.nodes[i]->SetFixed(m.p[i].x() < 1e-9);
        if (m.p[i].x() > L - 1e-9) {
            bool corner = m.p[i].y() < 1e-9 || m.p[i].y() > b - 1e-9;
            F[i].z() = P / ny * (corner ? 0.5 : 1.0);
        }
    }
    ChVector3d R = model.SolveStatic(F);
    std::cout << "cantilever strip: w/w_exact = " << model.MaxDeflection() / exact << std::endl;
    ASSERT_NEAR(model.MaxDeflection() / exact, 1.0, 0.02);
    ASSERT_NEAR(R.z(), -P, 1e-6 * P);
}
