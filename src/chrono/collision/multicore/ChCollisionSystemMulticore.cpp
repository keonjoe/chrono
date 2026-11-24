// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2021 projectchrono.org
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
// Chrono custom multicore collision system.
// Contains both the broadphase and the narrowphase methods.
//
// =============================================================================

#include <algorithm>

#include "chrono/physics/ChSystem.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChParticleCloud.h"
#include "chrono/utils/ChProfiler.h"

#include "chrono/multicore_math/thrust.h"

#include "chrono/collision/multicore/ChCollisionSystemMulticore.h"
#include "chrono/collision/multicore/ChRayTest.h"

namespace chrono {

CH_FACTORY_REGISTER(ChCollisionSystemMulticore)
CH_UPCASTING(ChCollisionSystemMulticore, ChCollisionSystem)

ChCollisionSystemMulticore::ChCollisionSystemMulticore() : use_aabb_active(false) {
    // Create the shared data structure with own state data
    cd_data = chrono_types::make_shared<ChCollisionData>(true);
    cd_data->collision_envelope = ChCollisionModel::GetDefaultSuggestedEnvelope();

    broadphase.cd_data = cd_data;
    narrowphase.cd_data = cd_data;
}

ChCollisionSystemMulticore::~ChCollisionSystemMulticore() {}

// -----------------------------------------------------------------------------

void ChCollisionSystemMulticore::SetEnvelope(double envelope) {
    cd_data->collision_envelope = real(envelope);
}

void ChCollisionSystemMulticore::SetBroadphaseGridResolution(const ChVector3i& num_bins) {
    broadphase.grid_resolution = vec3(num_bins.x(), num_bins.y(), num_bins.z());
    broadphase.grid_type = ChBroadphase::GridType::FIXED_RESOLUTION;
}

void ChCollisionSystemMulticore::SetBroadphaseGridSize(const ChVector3d& bin_size) {
    broadphase.bin_size = real3(bin_size.x(), bin_size.y(), bin_size.z());
    broadphase.grid_type = ChBroadphase::GridType::FIXED_RESOLUTION;
}

void ChCollisionSystemMulticore::SetBroadphaseGridDensity(double density) {
    broadphase.grid_density = real(density);
    broadphase.grid_type = ChBroadphase::GridType::FIXED_DENSITY;
}

void ChCollisionSystemMulticore::SetNarrowphaseAlgorithm(ChNarrowphase::Algorithm algorithm) {
    narrowphase.algorithm = algorithm;
}

void ChCollisionSystemMulticore::EnableActiveBoundingBox(const ChVector3d& aabb_min, const ChVector3d& aabb_max) {
    active_aabb_min = FromChVector(aabb_min);
    active_aabb_max = FromChVector(aabb_max);

    use_aabb_active = true;
}

void ChCollisionSystemMulticore::SetNumThreads(int nthreads) {
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif
}

// -----------------------------------------------------------------------------

void ChCollisionSystemMulticore::Add(std::shared_ptr<ChCollisionModel> model) {
    assert(!model->HasImplementation());

    auto ct_model = chrono_types::make_shared<ChCollisionModelMulticore>(model.get());
    ct_model->Populate();

    auto body_id = ct_model->GetBody()->GetIndex();
    short2 fam = S2(ct_model->model->GetFamilyGroup(), ct_model->model->GetFamilyMask());

    // The offset for this shape will the current total number of points in the convex data list
    auto& shape_data = cd_data->shape_data;
    int convex_data_offset = (int)shape_data.convex_rigid.size();

    // Insert the points into the global convex list
    shape_data.convex_rigid.insert(shape_data.convex_rigid.end(), ct_model->local_convex_data.begin(),
                                   ct_model->local_convex_data.end());

    // Shape index in the collision model
    int local_shape_index = 0;

    // Traverse all collision shapes in the model
    auto num_shapes = ct_model->m_shapes.size();
    assert(num_shapes == ct_model->m_ct_shapes.size());

    for (size_t i = 0; i < num_shapes; i++) {
        auto type = ct_model->m_shapes[i]->GetType();
        const auto& ct_shape = ct_model->m_ct_shapes[i];
        real3 obA = ct_shape->A;
        real3 obB = ct_shape->B;
        real3 obC = ct_shape->C;

        // Compute the global offset of the convex data structure based on the number of points already present
        int length = 1;
        int start;

        switch (type) {
            case ChCollisionShape::Type::SPHERE:
                start = (int)shape_data.sphere_rigid.size();
                shape_data.sphere_rigid.push_back(obB.x);
                break;
            case ChCollisionShape::Type::ELLIPSOID:
                start = (int)shape_data.box_like_rigid.size();
                shape_data.box_like_rigid.push_back(obB);
                break;
            case ChCollisionShape::Type::BOX:
                start = (int)shape_data.box_like_rigid.size();
                shape_data.box_like_rigid.push_back(obB);
                break;
            case ChCollisionShape::Type::CYLINDER:
                start = (int)shape_data.box_like_rigid.size();
                shape_data.box_like_rigid.push_back(obB);
                break;
            case ChCollisionShape::Type::CONE:
                start = (int)shape_data.box_like_rigid.size();
                shape_data.box_like_rigid.push_back(obB);
                break;
            case ChCollisionShape::Type::CAPSULE:
                start = (int)shape_data.capsule_rigid.size();
                shape_data.capsule_rigid.push_back(real2(obB.x, obB.y));
                break;
            case ChCollisionShape::Type::ROUNDEDBOX:
                start = (int)shape_data.rbox_like_rigid.size();
                shape_data.rbox_like_rigid.push_back(real4(obB, obC.x));
                break;
            case ChCollisionShape::Type::ROUNDEDCYL:
                start = (int)shape_data.rbox_like_rigid.size();
                shape_data.rbox_like_rigid.push_back(real4(obB, obC.x));
                break;
            case ChCollisionShape::Type::CONVEXHULL:
                start = (int)(obB.y + convex_data_offset);
                length = (int)obB.x;
                break;
            case ChCollisionShape::Type::TRIANGLE:
                start = (int)shape_data.triangle_rigid.size();
                shape_data.triangle_rigid.push_back(obA);
                shape_data.triangle_rigid.push_back(obB);
                shape_data.triangle_rigid.push_back(obC);
                break;
            default:
                start = -1;
                break;
        }

        shape_data.ObA_rigid.push_back(obA);
        shape_data.ObR_rigid.push_back(ct_shape->R);
        shape_data.start_rigid.push_back(start);
        shape_data.length_rigid.push_back(length);

        shape_data.fam_rigid.push_back(fam);
        shape_data.typ_rigid.push_back(type);
        shape_data.id_rigid.push_back(body_id);
        shape_data.local_rigid.push_back(local_shape_index);
        cd_data->num_rigid_shapes++;
        local_shape_index++;
    }

    ct_models.push_back(ct_model);
}

void ChCollisionSystemMulticore::Clear() {
    // Release all collision model implementations (their destructor detaches them from their models) and empty the
    // shape data.
    ct_models.clear();

    auto& sd = cd_data->shape_data;
    sd.fam_rigid.clear();
    sd.id_rigid.clear();
    sd.typ_rigid.clear();
    sd.local_rigid.clear();
    sd.start_rigid.clear();
    sd.length_rigid.clear();
    sd.ObR_rigid.clear();
    sd.ObA_rigid.clear();
    sd.obj_data_A_global.clear();
    sd.obj_data_R_global.clear();
    sd.sphere_rigid.clear();
    sd.box_like_rigid.clear();
    sd.triangle_rigid.clear();
    sd.capsule_rigid.clear();
    sd.rbox_like_rigid.clear();
    sd.convex_rigid.clear();
    sd.triangle_global.clear();
    cd_data->num_rigid_shapes = 0;
    cd_data->num_rigid_contacts = 0;
}

// Erase `count` entries starting at `start` from the given per-type shape dimension container and shift the start
// indices of the remaining shapes of the given `types` stored in that container.
template <typename T>
static void EraseShapeData(std::vector<T>& container,
                           int start,
                           int count,
                           shape_container& sd,
                           const std::vector<int>& types) {
    if (start < 0 || count <= 0 || start + count > (int)container.size())
        return;
    container.erase(container.begin() + start, container.begin() + start + count);
    for (size_t i = 0; i < sd.typ_rigid.size(); i++) {
        if (std::find(types.begin(), types.end(), sd.typ_rigid[i]) != types.end() && sd.start_rigid[i] > start)
            sd.start_rigid[i] -= count;
    }
}

void ChCollisionSystemMulticore::Remove(std::shared_ptr<ChCollisionModel> model) {
    if (!model || !model->HasImplementation())
        return;

    // Locate the multicore model and the range of its shapes in the shape data (models are stored in order of
    // insertion, each with a contiguous block of shapes).
    int shape_offset = 0;
    auto it = ct_models.begin();
    for (; it != ct_models.end(); ++it) {
        if ((*it)->model == model.get())
            break;
        shape_offset += (int)(*it)->m_shapes.size();
    }
    if (it == ct_models.end())
        return;  // Not in this collision system

    auto ct_model = *it;
    auto& sd = cd_data->shape_data;
    int num_shapes = (int)ct_model->m_shapes.size();
    if (shape_offset + num_shapes > (int)sd.typ_rigid.size())
        throw std::runtime_error("ChCollisionSystemMulticore::Remove: inconsistent collision shape data.");

    using T = ChCollisionShape::Type;
    static const std::vector<int> t_sphere = {(int)T::SPHERE};
    static const std::vector<int> t_box = {(int)T::ELLIPSOID, (int)T::BOX, (int)T::CYLINDER, (int)T::CONE};
    static const std::vector<int> t_capsule = {(int)T::CAPSULE};
    static const std::vector<int> t_rbox = {(int)T::ROUNDEDBOX, (int)T::ROUNDEDCYL};
    static const std::vector<int> t_triangle = {(int)T::TRIANGLE};
    static const std::vector<int> t_convex = {(int)T::CONVEXHULL};

    // Start of the block of convex hull points of this model (inserted contiguously by Add).
    int convex_start = -1;
    for (int k = 0; k < num_shapes; k++) {
        int i = shape_offset + k;
        if (sd.typ_rigid[i] == (int)T::CONVEXHULL) {
            convex_start = sd.start_rigid[i] - (int)ct_model->m_ct_shapes[k]->B.y;
            break;
        }
    }

    // Remove the per-type dimension data of the model shapes (last first, so that start indices stay valid).
    // The shapes of the model are first flagged with an invalid type so that the start index shift only applies to
    // the shapes that remain.
    std::vector<std::pair<int, int>> removed;  // (type, start)
    for (int k = 0; k < num_shapes; k++) {
        int i = shape_offset + k;
        removed.push_back({sd.typ_rigid[i], sd.start_rigid[i]});
        sd.typ_rigid[i] = -1;
    }
    std::sort(removed.begin(), removed.end(),
              [](const std::pair<int, int>& a, const std::pair<int, int>& b) { return a.second > b.second; });
    for (const auto& r : removed) {
        int type = r.first;
        int start = r.second;
        if (type == (int)T::SPHERE)
            EraseShapeData(sd.sphere_rigid, start, 1, sd, t_sphere);
        else if (std::find(t_box.begin(), t_box.end(), type) != t_box.end())
            EraseShapeData(sd.box_like_rigid, start, 1, sd, t_box);
        else if (type == (int)T::CAPSULE)
            EraseShapeData(sd.capsule_rigid, start, 1, sd, t_capsule);
        else if (std::find(t_rbox.begin(), t_rbox.end(), type) != t_rbox.end())
            EraseShapeData(sd.rbox_like_rigid, start, 1, sd, t_rbox);
        else if (type == (int)T::TRIANGLE)
            EraseShapeData(sd.triangle_rigid, start, 3, sd, t_triangle);
    }
    if (convex_start >= 0)
        EraseShapeData(sd.convex_rigid, convex_start, (int)ct_model->local_convex_data.size(), sd, t_convex);

    // Remove the per-shape data.
    auto erase_range = [shape_offset, num_shapes](auto& v) {
        if ((int)v.size() >= shape_offset + num_shapes)
            v.erase(v.begin() + shape_offset, v.begin() + shape_offset + num_shapes);
    };
    erase_range(sd.fam_rigid);
    erase_range(sd.id_rigid);
    erase_range(sd.typ_rigid);
    erase_range(sd.local_rigid);
    erase_range(sd.start_rigid);
    erase_range(sd.length_rigid);
    erase_range(sd.ObR_rigid);
    erase_range(sd.ObA_rigid);
    erase_range(sd.obj_data_A_global);
    erase_range(sd.obj_data_R_global);
    sd.triangle_global.resize(sd.triangle_rigid.size());
    cd_data->num_rigid_shapes -= num_shapes;

    // Contacts found in the last collision detection may refer to removed shapes.
    cd_data->num_rigid_contacts = 0;

    ct_models.erase(it);
    model->RemoveImplementation();
}

// Body indices change when bodies are removed from the system: refresh the body ID of every shape.
void ChCollisionSystemMulticore::UpdateShapeBodyIDs() {
    auto& id_rigid = cd_data->shape_data.id_rigid;
    size_t i = 0;
    for (const auto& ct_model : ct_models) {
        uint body_id = (uint)ct_model->GetBody()->GetIndex();
        for (size_t k = 0; k < ct_model->m_shapes.size() && i < id_rigid.size(); k++, i++)
            id_rigid[i] = body_id;
    }
}

// -----------------------------------------------------------------------------

bool ChCollisionSystemMulticore::GetActiveBoundingBox(ChVector3d& aabb_min, ChVector3d& aabb_max) const {
    aabb_min = ToChVector(active_aabb_min);
    aabb_max = ToChVector(active_aabb_max);

    return use_aabb_active;
}

ChAABB ChCollisionSystemMulticore::GetBoundingBox() const {
    ChVector3d aabb_min((double)cd_data->min_bounding_point.x, (double)cd_data->min_bounding_point.y,
                        (double)cd_data->min_bounding_point.z);
    ChVector3d aabb_max((double)cd_data->max_bounding_point.x, (double)cd_data->max_bounding_point.y,
                        (double)cd_data->max_bounding_point.z);

    return ChAABB(aabb_min, aabb_max);
}

void ChCollisionSystemMulticore::ResetTimers() {
    m_timer_broad.reset();
    m_timer_narrow.reset();
}

double ChCollisionSystemMulticore::GetTimerCollisionBroad() const {
    return m_timer_broad();
}

double ChCollisionSystemMulticore::GetTimerCollisionNarrow() const {
    return m_timer_narrow();
}

// -----------------------------------------------------------------------------

void ChCollisionSystemMulticore::PreProcess() {
    assert(cd_data->owns_state_data);

    std::vector<real3>& position = *cd_data->state_data.pos_rigid;
    std::vector<quaternion>& rotation = *cd_data->state_data.rot_rigid;
    std::vector<char>& active = *cd_data->state_data.active_rigid;
    std::vector<char>& collide = *cd_data->state_data.collide_rigid;

    const auto& blist = m_system->GetBodies();
    int m_num_bodies = static_cast<int>(blist.size());

    position.resize(m_num_bodies);
    rotation.resize(m_num_bodies);
    active.resize(m_num_bodies);
    collide.resize(m_num_bodies);

    cd_data->state_data.num_rigid_bodies = m_num_bodies;
    cd_data->state_data.num_particles = 0;

    UpdateShapeBodyIDs();

#pragma omp parallel for
    for (int i = 0; i < m_num_bodies; i++) {
        const auto& body = blist[i];

        const ChVector3d& body_pos = body->GetPos();
        const ChQuaternion<>& body_rot = body->GetRot();

        position[i] = real3(body_pos.x(), body_pos.y(), body_pos.z());
        rotation[i] = quaternion(body_rot.e0(), body_rot.e1(), body_rot.e2(), body_rot.e3());

        active[i] = body->IsActive();
        collide[i] = body->IsCollisionEnabled();
    }
}

void ChCollisionSystemMulticore::PostProcess() {
    if (use_aabb_active) {
        const auto& active = *cd_data->state_data.active_rigid;
        auto& blist = m_system->GetBodies();
        int m_num_bodies = static_cast<int>(blist.size());

#pragma omp parallel for
        for (int i = 0; i < m_num_bodies; i++) {
            blist[i]->SetSleeping(!active[i]);
        }
    }
}

void ChCollisionSystemMulticore::Run() {
    ResetTimers();

    if (use_aabb_active) {
        std::vector<char>& active = *cd_data->state_data.active_rigid;
        const std::vector<char>& collide = *cd_data->state_data.collide_rigid;

        body_active.resize(cd_data->state_data.num_rigid_bodies);
        std::fill(body_active.begin(), body_active.end(), false);

        GetOverlappingAABB(body_active, active_aabb_min, active_aabb_max);

#pragma omp parallel for
        for (int i = 0; i < active.size(); i++) {
            if (active[i] != 0 && collide[i] != 0) {
                active[i] = body_active[i];
            }
        }
    }

    // Broadphase
    {
        CH_PROFILE("Broad-phase");
        m_timer_broad.start();
        GenerateAABB();
        broadphase.Process();
        m_timer_broad.stop();
    }

    // Narrowphase
    {
        CH_PROFILE("Narrow-phase");
        m_timer_narrow.start();
        narrowphase.Process();
        m_timer_narrow.stop();
    }
}

// -----------------------------------------------------------------------------

void ChCollisionSystemMulticore::ReportContacts(ChContactContainer* contact_container) {
    const auto& blist = m_system->GetBodies();

    // Resize global arrays with composite material properties.
    // NOTE: important to do this here, to set size to zero if no contacts (in case some other added by a custom user
    // callback)
    contact_container->BeginAddContact();

    const auto& bids = cd_data->bids_rigid_rigid;          // global IDs of bodies in contact
    const auto& sids = cd_data->contact_shapeIDs;          // global IDs of shapes in contact
    const auto& sindex = cd_data->shape_data.local_rigid;  // collision model indexes of shapes in contact

    // Loop over all current contacts, create the cinfo structure and add contact to the container.
    // Note that inclusions in the contact container cannot be done in parallel.
    for (uint i = 0; i < cd_data->num_rigid_contacts; i++) {
        auto b1 = bids[i].x;                  // global IDs of bodies in contact
        auto b2 = bids[i].y;                  //
        auto s1 = int(sids[i] >> 32);         // global IDs of shapes in contact
        auto s2 = int(sids[i] & 0xffffffff);  //
        auto s1_index = sindex[s1];           // indexes of shapes in contact within their collision model
        auto s2_index = sindex[s2];           //

        auto modelA = (ChCollisionModelMulticore*)blist[b1]->GetCollisionModel()->GetImplementation();
        auto modelB = (ChCollisionModelMulticore*)blist[b2]->GetCollisionModel()->GetImplementation();

        ChCollisionInfo cinfo;
        cinfo.modelA = blist[b1]->GetCollisionModel().get();
        cinfo.modelB = blist[b2]->GetCollisionModel().get();
        cinfo.shapeA = modelA->m_shapes[s1_index].get();
        cinfo.shapeB = modelB->m_shapes[s2_index].get();
        cinfo.vN = ToChVector(cd_data->norm_rigid_rigid[i]);
        cinfo.vpA = ToChVector(cd_data->cpta_rigid_rigid[i]);
        cinfo.vpB = ToChVector(cd_data->cptb_rigid_rigid[i]);
        cinfo.distance = cd_data->dpth_rigid_rigid[i];
        cinfo.eff_radius = cd_data->erad_rigid_rigid[i];

        // Execute user custom callback, if any
        bool add_contact = true;
        if (this->narrow_callback)
            add_contact = this->narrow_callback->OnNarrowphase(cinfo);

        if (add_contact)
            contact_container->AddContact(cinfo);
    }

    contact_container->EndAddContact();
}

// -----------------------------------------------------------------------------

static void ComputeAABBSphere(const real& radius,
                              const real3& lpositon,
                              const real3& position,
                              const quaternion& body_rotation,
                              real3& minp,
                              real3& maxp) {
    real3 pos = Rotate(lpositon, body_rotation) + position;
    minp = pos - radius;
    maxp = pos + radius;
}

static void ComputeAABBTriangle(const real3& A, const real3& B, const real3& C, real3& minp, real3& maxp) {
    minp.x = std::min(A.x, std::min(B.x, C.x));
    minp.y = std::min(A.y, std::min(B.y, C.y));
    minp.z = std::min(A.z, std::min(B.z, C.z));
    maxp.x = std::max(A.x, std::max(B.x, C.x));
    maxp.y = std::max(A.y, std::max(B.y, C.y));
    maxp.z = std::max(A.z, std::max(B.z, C.z));
}

static void ComputeAABBBox(const real3& dim,
                           const real3& lposition,
                           const real3& position,
                           const quaternion& rotation,
                           const quaternion& body_rotation,
                           real3& minp,
                           real3& maxp) {
    real3 temp = AbsRotate(rotation, dim);
    real3 pos = Rotate(lposition, body_rotation) + position;
    minp = pos - temp;
    maxp = pos + temp;
}

/*
static void ComputeAABBCone(const real3& dim,
                            const real3& lposition,
                            const real3& position,
                            const quaternion& rotation,
                            const quaternion& body_rotation,
                            real3& minp,
                            real3& maxp) {
    real3 temp = AbsRotate(rotation, real3(dim.x, dim.y, dim.z / 2.0));
    real3 pos = Rotate(lposition - real3(0, 0, dim.z / 2.0), body_rotation) + position;
    minp = pos - temp;
    maxp = pos + temp;
}
*/

static void ComputeAABBConvex(const real3* convex_points,
                              const int start,
                              const int size,
                              const real3& lpos,
                              const real3& pos,
                              const quaternion& rot,
                              real3& minp,
                              real3& maxp) {
    real3 point_0 = Rotate(convex_points[start] + lpos, rot) + pos;

    minp = maxp = point_0;
    for (int i = start; i < start + size; i++) {
        real3 p = Rotate(convex_points[i] + lpos, rot) + pos;
        if (minp.x > p.x) {
            minp.x = p.x;
        }
        if (minp.y > p.y) {
            minp.y = p.y;
        }
        if (minp.z > p.z) {
            minp.z = p.z;
        }
        if (maxp.x < p.x) {
            maxp.x = p.x;
        }
        if (maxp.y < p.y) {
            maxp.y = p.y;
        }
        if (maxp.z < p.z) {
            maxp.z = p.z;
        }
    }
}

void ChCollisionSystemMulticore::GenerateAABB() {
    if (cd_data->num_rigid_shapes > 0) {
        const real envelope = cd_data->collision_envelope;
        const std::vector<shape_type>& typ_rigid = cd_data->shape_data.typ_rigid;
        const std::vector<int>& start_rigid = cd_data->shape_data.start_rigid;
        const std::vector<uint>& id_rigid = cd_data->shape_data.id_rigid;
        const std::vector<real3>& obj_data_A = cd_data->shape_data.ObA_rigid;
        const std::vector<quaternion>& obj_data_R = cd_data->shape_data.ObR_rigid;
        const std::vector<real3>& convex_rigid = cd_data->shape_data.convex_rigid;

        const std::vector<real3>& pos_rigid = *cd_data->state_data.pos_rigid;
        const std::vector<quaternion>& body_rot = *cd_data->state_data.rot_rigid;

        const uint num_rigid_shapes = cd_data->num_rigid_shapes;

        std::vector<real3>& aabb_min = cd_data->aabb_min;
        std::vector<real3>& aabb_max = cd_data->aabb_max;

        aabb_min.resize(num_rigid_shapes);
        aabb_max.resize(num_rigid_shapes);

#pragma omp parallel for
        for (int index = 0; index < (signed)num_rigid_shapes; index++) {
            // Shape data
            shape_type type = typ_rigid[index];
            real3 local_pos = obj_data_A[index];
            quaternion local_rot = obj_data_R[index];
            uint id = id_rigid[index];  // The rigid body corresponding to this shape
            int start = start_rigid[index];

            // Body data
            if (id == UINT_MAX)
                continue;

            real3 position = pos_rigid[id];
            quaternion rotation = Mult(body_rot[id], local_rot);
            real3 temp_min;
            real3 temp_max;

            if (type == ChCollisionShape::Type::SPHERE) {
                real radius = cd_data->shape_data.sphere_rigid[start];
                ComputeAABBSphere(radius + envelope, local_pos, position, body_rot[id], temp_min, temp_max);

            } else if (type == ChCollisionShape::Type::ELLIPSOID || type == ChCollisionShape::Type::BOX ||
                       type == ChCollisionShape::Type::CYLINDER || type == ChCollisionShape::Type::CONE) {
                real3 B = cd_data->shape_data.box_like_rigid[start];
                ComputeAABBBox(B + envelope, local_pos, position, rotation, body_rot[id], temp_min, temp_max);

            } else if (type == ChCollisionShape::Type::ROUNDEDBOX || type == ChCollisionShape::Type::ROUNDEDCYL) {
                real4 T = cd_data->shape_data.rbox_like_rigid[start];
                real3 B = real3(T.x, T.y, T.z) + envelope;
                ComputeAABBBox(B, local_pos, position, rotation, body_rot[id], temp_min, temp_max);

            } else if (type == ChCollisionShape::Type::CAPSULE) {
                real2 T = cd_data->shape_data.capsule_rigid[start];
                real3 B = real3(T.x, T.x + T.y, T.x) + envelope;
                ComputeAABBBox(B, local_pos, position, rotation, body_rot[id], temp_min, temp_max);

            } else if (type == ChCollisionShape::Type::CONVEXHULL) {
                int length = cd_data->shape_data.length_rigid[index];
                ComputeAABBConvex(convex_rigid.data(), start, length, local_pos, position, rotation, temp_min,
                                  temp_max);
                temp_min -= envelope;
                temp_max += envelope;

            } else if (type == ChCollisionShape::Type::TRIANGLE) {
                real3 A, B, C;

                A = cd_data->shape_data.triangle_rigid[start + 0];
                B = cd_data->shape_data.triangle_rigid[start + 1];
                C = cd_data->shape_data.triangle_rigid[start + 2];

                A = Rotate(A, body_rot[id]) + position;
                B = Rotate(B, body_rot[id]) + position;
                C = Rotate(C, body_rot[id]) + position;

                ComputeAABBTriangle(A, B, C, temp_min, temp_max);

            } else {
                continue;
            }

            aabb_min[index] = temp_min;
            aabb_max[index] = temp_max;
        }
    }
}

void ChCollisionSystemMulticore::GetOverlappingAABB(std::vector<char>& active_id, real3 Amin, real3 Amax) {
    GenerateAABB();

#pragma omp parallel for
    for (int i = 0; i < cd_data->shape_data.typ_rigid.size(); i++) {
        real3 Bmin = cd_data->aabb_min[i];
        real3 Bmax = cd_data->aabb_max[i];

        bool inContact = (Amin.x <= Bmax.x && Bmin.x <= Amax.x) && (Amin.y <= Bmax.y && Bmin.y <= Amax.y) &&
                         (Amin.z <= Bmax.z && Bmin.z <= Amax.z);
        if (inContact) {
            active_id[cd_data->shape_data.id_rigid[i]] = true;
        }
    }
}

std::vector<vec2> ChCollisionSystemMulticore::GetOverlappingPairs() {
    std::vector<vec2> pairs;
    pairs.resize(cd_data->pair_shapeIDs.size());
    for (int i = 0; i < cd_data->pair_shapeIDs.size(); i++) {
        vec2 pair = I2(int(cd_data->pair_shapeIDs[i] >> 32), int(cd_data->pair_shapeIDs[i] & 0xffffffff));
        pairs[i] = pair;
    }
    return pairs;
}

// -----------------------------------------------------------------------------

bool ChCollisionSystemMulticore::RayHit(const ChVector3d& from, const ChVector3d& to, ChRayhitResult& result) const {
    if (cd_data->num_active_bins == 0) {
        result.hit = false;
        return false;
    }

    ChRayTest tester(cd_data);
    ChRayTest::RayHitInfo info;
    if (tester.Check(FromChVector(from), FromChVector(to), info)) {
        // Hit point
        result.hit = true;
        result.abs_hitNormal = ToChVector(info.normal);
        result.abs_hitPoint = ToChVector(info.point);
        result.dist_factor = info.t;

        // ID of the body carring the closest hit shape
        uint bid = cd_data->shape_data.id_rigid[info.shapeID];

        // Collision model of hit body
        result.hitModel = m_system->GetBodies()[bid]->GetCollisionModel().get();

        return true;
    }

    result.hit = false;
    return false;
}

bool ChCollisionSystemMulticore::RayHit(const ChVector3d& from,
                                        const ChVector3d& to,
                                        ChCollisionModel* model,
                                        ChRayhitResult& result) const {
    result.hit = false;
    if (!model)
        return false;

    // Shapes of the given collision model.
    std::vector<int> shapes;
    int shape_offset = 0;
    for (const auto& ct_model : ct_models) {
        int n = (int)ct_model->m_shapes.size();
        if (ct_model->model == model) {
            for (int k = 0; k < n; k++)
                shapes.push_back(shape_offset + k);
            break;
        }
        shape_offset += n;
    }
    if (shapes.empty())
        return false;

    ChRayTest tester(cd_data);
    ChRayTest::RayHitInfo info;
    if (!tester.CheckShapes(FromChVector(from), FromChVector(to), shapes, info))
        return false;

    result.hit = true;
    result.abs_hitNormal = ToChVector(info.normal);
    result.abs_hitPoint = ToChVector(info.point);
    result.dist_factor = info.t;
    result.hitModel = model;
    return true;
}

// -----------------------------------------------------------------------------

void DrawHemisphere(ChCollisionSystem::VisualizationCallback* vis,
                    const ChCoordsys<>& csys,
                    double radius,
                    const ChColor& color) {
    int vstep = 30;  // degrees
    int rstep = 30;  // degrees

    for (int j = 0; j < 90; j += vstep) {
        double y_low = radius * std::sin(j * CH_DEG_TO_RAD);
        double r_low = radius * std::cos(j * CH_DEG_TO_RAD);

        double y_top = radius * std::sin((j + vstep) * CH_DEG_TO_RAD);
        double r_top = radius * std::cos((j + vstep) * CH_DEG_TO_RAD);

        ChVector3d crt(r_low, y_low, 0);
        for (int i = 0; i < 360; i += rstep) {
            double rc = r_top * std::cos(i * CH_DEG_TO_RAD);
            double rs = r_top * std::sin(i * CH_DEG_TO_RAD);
            ChVector3d up(rc, y_top, rs);
            vis->DrawLine(csys.TransformPointLocalToParent(crt), csys.TransformPointLocalToParent(up), color);
            rc = r_low * std::cos((i + rstep) * CH_DEG_TO_RAD);
            rs = r_low * std::sin((i + rstep) * CH_DEG_TO_RAD);
            ChVector3d next(rc, y_low, rs);
            vis->DrawLine(csys.TransformPointLocalToParent(crt), csys.TransformPointLocalToParent(next), color);
            crt = next;
        }
    }
}

void DrawSphere(ChCollisionSystem::VisualizationCallback* vis,
                const ChCoordsys<>& csys,
                double radius,
                const ChColor& color) {
    DrawHemisphere(vis, csys, radius, color);
    ChCoordsys<> csys1 = csys;
    csys1.rot = csys.rot * QuatFromAngleX(CH_PI);
    DrawHemisphere(vis, csys1, radius, color);
}

void DrawBox(ChCollisionSystem::VisualizationCallback* vis,
             const ChCoordsys<>& csys,
             const ChVector3d& hdim,
             const ChColor& color) {
    vis->DrawLine(csys.TransformPointLocalToParent(ChVector3d(-hdim[0], -hdim[1], -hdim[2])),
                  csys.TransformPointLocalToParent(ChVector3d(+hdim[0], -hdim[1], -hdim[2])), color);
    vis->DrawLine(csys.TransformPointLocalToParent(ChVector3d(+hdim[0], -hdim[1], -hdim[2])),
                  csys.TransformPointLocalToParent(ChVector3d(+hdim[0], +hdim[1], -hdim[2])), color);
    vis->DrawLine(csys.TransformPointLocalToParent(ChVector3d(+hdim[0], +hdim[1], -hdim[2])),
                  csys.TransformPointLocalToParent(ChVector3d(-hdim[0], +hdim[1], -hdim[2])), color);
    vis->DrawLine(csys.TransformPointLocalToParent(ChVector3d(-hdim[0], +hdim[1], -hdim[2])),
                  csys.TransformPointLocalToParent(ChVector3d(-hdim[0], -hdim[1], -hdim[2])), color);
    vis->DrawLine(csys.TransformPointLocalToParent(ChVector3d(-hdim[0], -hdim[1], -hdim[2])),
                  csys.TransformPointLocalToParent(ChVector3d(-hdim[0], -hdim[1], +hdim[2])), color);
    vis->DrawLine(csys.TransformPointLocalToParent(ChVector3d(+hdim[0], -hdim[1], -hdim[2])),
                  csys.TransformPointLocalToParent(ChVector3d(+hdim[0], -hdim[1], +hdim[2])), color);
    vis->DrawLine(csys.TransformPointLocalToParent(ChVector3d(+hdim[0], +hdim[1], -hdim[2])),
                  csys.TransformPointLocalToParent(ChVector3d(+hdim[0], +hdim[1], +hdim[2])), color);
    vis->DrawLine(csys.TransformPointLocalToParent(ChVector3d(-hdim[0], +hdim[1], -hdim[2])),
                  csys.TransformPointLocalToParent(ChVector3d(-hdim[0], +hdim[1], +hdim[2])), color);
    vis->DrawLine(csys.TransformPointLocalToParent(ChVector3d(-hdim[0], -hdim[1], +hdim[2])),
                  csys.TransformPointLocalToParent(ChVector3d(+hdim[0], -hdim[1], +hdim[2])), color);
    vis->DrawLine(csys.TransformPointLocalToParent(ChVector3d(+hdim[0], -hdim[1], +hdim[2])),
                  csys.TransformPointLocalToParent(ChVector3d(+hdim[0], +hdim[1], +hdim[2])), color);
    vis->DrawLine(csys.TransformPointLocalToParent(ChVector3d(+hdim[0], +hdim[1], +hdim[2])),
                  csys.TransformPointLocalToParent(ChVector3d(-hdim[0], +hdim[1], +hdim[2])), color);
    vis->DrawLine(csys.TransformPointLocalToParent(ChVector3d(-hdim[0], +hdim[1], +hdim[2])),
                  csys.TransformPointLocalToParent(ChVector3d(-hdim[0], -hdim[1], +hdim[2])), color);
}

void DrawCylinder(ChCollisionSystem::VisualizationCallback* vis,
                  const ChCoordsys<>& csys,
                  double radius,
                  double hlen,
                  const ChColor& color) {
    int rstep = 30;
    ChVector3d p_start(radius, 0, -hlen);
    ChVector3d p_end(radius, 0, +hlen);

    for (int i = 0; i < 360; i += rstep) {
        vis->DrawLine(csys.TransformPointLocalToParent(p_start), csys.TransformPointLocalToParent(p_end), color);
        double rc = radius * std::cos((i + rstep) * CH_DEG_TO_RAD);
        double rs = radius * std::sin((i + rstep) * CH_DEG_TO_RAD);
        ChVector3d start(rc, rs, -hlen);
        ChVector3d end(rc, rs, +hlen);
        vis->DrawLine(csys.TransformPointLocalToParent(p_start), csys.TransformPointLocalToParent(start), color);
        vis->DrawLine(csys.TransformPointLocalToParent(p_end), csys.TransformPointLocalToParent(end), color);
        p_start = start;
        p_end = end;
    }
}

void DrawCone(ChCollisionSystem::VisualizationCallback* vis,
              const ChCoordsys<>& csys,
              double radius,
              double hlen,
              const ChColor& color) {
    int rstep = 30;
    ChVector3d p_start(radius, 0, 0);
    ChVector3d p_end(0, 0, 2 * hlen);

    for (int i = 0; i < 360; i += rstep) {
        vis->DrawLine(csys.TransformPointLocalToParent(p_start), csys.TransformPointLocalToParent(p_end), color);
        double rc = radius * std::cos((i + rstep) * CH_DEG_TO_RAD);
        double rs = radius * std::sin((i + rstep) * CH_DEG_TO_RAD);
        ChVector3d start(rc, rs, 0);
        ChVector3d end(rc, rs, 2 * hlen);
        vis->DrawLine(csys.TransformPointLocalToParent(p_start), csys.TransformPointLocalToParent(start), color);
        p_start = start;
    }
}

void DrawCapsule(ChCollisionSystem::VisualizationCallback* vis,
                 const ChCoordsys<>& csys,
                 double radius,
                 double hlen,
                 const ChColor& color) {
    int rstep = 30;

    for (int i = 0; i < 360; i += rstep) {
        double rc = radius * std::cos((i + rstep) * CH_DEG_TO_RAD);
        double rs = radius * std::sin((i + rstep) * CH_DEG_TO_RAD);
        ChVector3d start(rc, rs, -hlen);
        ChVector3d end(rc, rs, +hlen);
        vis->DrawLine(csys.TransformPointLocalToParent(start), csys.TransformPointLocalToParent(end), color);
    }
    ChCoordsys<> csys1;
    csys1.pos = csys.TransformPointLocalToParent(ChVector3d(0, 0, hlen));
    csys1.rot = csys.rot;
    DrawHemisphere(vis, csys1, radius, color);
    ChCoordsys<> csys2;
    csys2.pos = csys.TransformPointLocalToParent(ChVector3d(0, 0, -hlen));
    csys2.rot = csys.rot * QuatFromAngleX(CH_PI);
    DrawHemisphere(vis, csys2, radius, color);
}

void ChCollisionSystemMulticore::Visualize(int flags) {
    if (!vis_callback || flags == VisualizationModes::VIS_None)
        return;

    if (flags & VisualizationModes::VIS_Shapes)
        VisualizeShapes();
    if (flags & VisualizationModes::VIS_Aabb)
        VisualizeAABB();
    if (flags & VisualizationModes::VIS_Contacts)
        VisualizeContacts();
}

void ChCollisionSystemMulticore::VisualizeShapes() {
    const real envelope = cd_data->collision_envelope;
    const std::vector<shape_type>& typ_rigid = cd_data->shape_data.typ_rigid;
    const std::vector<int>& start_rigid = cd_data->shape_data.start_rigid;
    const std::vector<uint>& id_rigid = cd_data->shape_data.id_rigid;
    const std::vector<real3>& obj_data_A = cd_data->shape_data.ObA_rigid;
    const std::vector<quaternion>& obj_data_R = cd_data->shape_data.ObR_rigid;

    const std::vector<real3>& pos_rigid = *cd_data->state_data.pos_rigid;
    const std::vector<quaternion>& body_rot = *cd_data->state_data.rot_rigid;

    const uint num_rigid_shapes = cd_data->num_rigid_shapes;

    for (uint index = 0; index < num_rigid_shapes; index++) {
        // Shape data
        shape_type type = typ_rigid[index];
        real3 local_pos = obj_data_A[index];
        quaternion local_rot = obj_data_R[index];
        uint id = id_rigid[index];  // The rigid body corresponding to this shape
        int start = start_rigid[index];

        // Body data
        if (id == UINT_MAX)
            continue;

        real3 position = pos_rigid[id] + Rotate(local_pos, body_rot[id]);
        quaternion rotation = Mult(body_rot[id], local_rot);

        switch (type) {
            case ChCollisionShape::Type::SPHERE: {
                real radius = cd_data->shape_data.sphere_rigid[start];
                DrawSphere(vis_callback.get(), ChCoordsys<>(ToChVector(position), ToChQuaternion(rotation)),
                           double(radius + envelope), ChColor(1, 0, 0));
                break;
            }
            case ChCollisionShape::Type::BOX: {
                const real3& B = cd_data->shape_data.box_like_rigid[start];
                DrawBox(vis_callback.get(), ChCoordsys<>(ToChVector(position), ToChQuaternion(rotation)),
                        ToChVector(B + envelope), ChColor(1, 0, 0));
                break;
            }
            case ChCollisionShape::Type::CYLINDER: {
                const real3& B = cd_data->shape_data.box_like_rigid[start];
                DrawCylinder(vis_callback.get(), ChCoordsys<>(ToChVector(position), ToChQuaternion(rotation)),
                             double(B.x + envelope), double(B.z + envelope), ChColor(1, 0, 0));
                break;
            }
            case ChCollisionShape::Type::CAPSULE: {
                const real3& B = cd_data->shape_data.box_like_rigid[start];
                DrawCapsule(vis_callback.get(), ChCoordsys<>(ToChVector(position), ToChQuaternion(rotation)),
                            double(B.x + envelope), double(B.y + envelope), ChColor(1, 0, 0));
                break;
            }
            case ChCollisionShape::Type::ELLIPSOID: {
                break;
            }
            case ChCollisionShape::Type::CONE: {
                const real3& B = cd_data->shape_data.box_like_rigid[start];
                DrawCone(vis_callback.get(), ChCoordsys<>(ToChVector(position), ToChQuaternion(rotation)),
                         double(B.x + envelope), double(B.y + envelope), ChColor(1, 0, 0));
                break;
            }
            case ChCollisionShape::Type::TRIANGLE: {
                real3 A = Rotate(cd_data->shape_data.triangle_rigid[start + 0], body_rot[id]) + pos_rigid[id];
                real3 B = Rotate(cd_data->shape_data.triangle_rigid[start + 1], body_rot[id]) + pos_rigid[id];
                real3 C = Rotate(cd_data->shape_data.triangle_rigid[start + 2], body_rot[id]) + pos_rigid[id];
                vis_callback->DrawLine(ToChVector(A), ToChVector(B), ChColor(1, 0, 0));
                vis_callback->DrawLine(ToChVector(B), ToChVector(C), ChColor(1, 0, 0));
                vis_callback->DrawLine(ToChVector(C), ToChVector(A), ChColor(1, 0, 0));
                break;
            }
        }
    }
}

void ChCollisionSystemMulticore::VisualizeAABB() {
    const uint num_rigid_shapes = cd_data->num_rigid_shapes;
    std::vector<real3>& aabb_min = cd_data->aabb_min;
    std::vector<real3>& aabb_max = cd_data->aabb_max;

    for (uint index = 0; index < num_rigid_shapes; index++) {
        real3 center = cd_data->global_origin + 0.5 * (aabb_max[index] + aabb_min[index]);
        real3 hdim = 0.5 * (aabb_max[index] - aabb_min[index]);
        DrawBox(vis_callback.get(), ChCoordsys<>(ToChVector(center), QUNIT), ToChVector(hdim), ChColor(0, 0, 1));
    }
}

void ChCollisionSystemMulticore::VisualizeContacts() {
    for (uint i = 0; i < cd_data->num_rigid_contacts; i++) {
        real3 from = cd_data->cptb_rigid_rigid[i];
        real3 to = cd_data->cptb_rigid_rigid[i] + vis_callback->GetNormalScale() * cd_data->norm_rigid_rigid[i];
        vis_callback->DrawLine(ToChVector(from), ToChVector(to), ChColor(1, 1, 0));
    }
}

void ChCollisionSystemMulticore::ArchiveOut(ChArchiveOut& archive_out) {
    // version number
    archive_out.VersionWrite<ChCollisionSystemMulticore>();
    // serialize parent class
    ChCollisionSystem::ArchiveOut(archive_out);
}

void ChCollisionSystemMulticore::ArchiveIn(ChArchiveIn& archive_in) {
    // version number
    /*int version =*/archive_in.VersionRead<ChCollisionSystemMulticore>();
    // deserialize parent class
    ChCollisionSystem::ArchiveIn(archive_in);
}

}  // end namespace chrono
