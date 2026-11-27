/*
***CHRONO***
*/

#include <cmath>

#include "cbtBarrelShape.h"

#include "BulletCollision/CollisionShapes/cbtCollisionMargin.h"
#include "LinearMath/cbtQuaternion.h"

cbtBarrelShape::cbtBarrelShape(cbtScalar sY_low, cbtScalar sY_high, cbtScalar sR_vert, cbtScalar sR_hor, cbtScalar sR_offset, cbtScalar offset) {
    Y_low = sY_low;
    Y_high = sY_high;
    R_vert = sR_vert;
    R_hor = sR_hor;
    R_offset = sR_offset;
    m_offset = offset;
    m_shapeType = BARREL_SHAPE_PROXYTYPE;
}

// Radius at level y of the lathed ellipse offset inward by s: the offset point of the ellipse point at angle t is
// (R_hor cos t, R_vert sin t) - s * normal(t); its height increases with t, so t is found by bisection.
cbtScalar cbtBarrelShape::getCoreRadius(cbtScalar y, cbtScalar s) const {
    cbtScalar t_lo = -SIMD_HALF_PI;
    cbtScalar t_hi = SIMD_HALF_PI;
    cbtScalar r = 0;
    for (int i = 0; i < 40; i++) {
        cbtScalar t = (t_lo + t_hi) / 2;
        cbtScalar nx = std::cos(t) / R_hor;
        cbtScalar ny = std::sin(t) / R_vert;
        cbtScalar nl = std::sqrt(nx * nx + ny * ny);
        cbtScalar ty = R_vert * std::sin(t) - s * ny / nl;
        r = R_hor * std::cos(t) - s * nx / nl;
        if (ty < y)
            t_lo = t;
        else
            t_hi = t;
    }
    return r + R_offset;
}

#include <stdio.h>
cbtVector3 cbtBarrelShape::localGetSupportingVertexWithoutMargin(const cbtVector3& vec0) const {
    cbtVector3 supVec(0, 0, 0);
    cbtVector3 supVecD;

    // The shape with its margin is the barrel inflated by the offset. Without the margin it is the barrel shrunk by
    // s = margin - offset: the lathed ellipse is offset inward by s (exact where its radius of curvature exceeds s)
    // and the discs are moved inward by s, so that the edges are rounded with radius s.
    cbtScalar s = getMargin() - m_offset;
    cbtScalar sY_high = Y_high - s;
    cbtScalar sY_low = Y_low + s;

    cbtVector3 dir = vec0;
    cbtScalar lenSqr = dir.length2();
    if (lenSqr < SIMD_EPSILON * SIMD_EPSILON)
        dir.setValue(1, 0, 0);
    else
        dir *= cbtScalar(1) / cbtSqrt(lenSqr);

    cbtScalar pY = dir.y();
    cbtScalar pR = std::sqrt(dir.z() * dir.z() + dir.x() * dir.x());

    // direction along the axis: center of the top or bottom disc (or of the top or bottom of the ellipse)
    if (pR < SIMD_EPSILON) {
        if (pY > 0)
            supVec.setValue(0, (std::abs(Y_high) < R_vert) ? sY_high : R_vert - s, 0);
        else
            supVec.setValue(0, (std::abs(Y_low) < R_vert) ? sY_low : -R_vert + s, 0);
        return supVec;
    }

    // support point on the lathed ellipse?
    cbtScalar pH = pR;
    cbtScalar dX = dir.x() / pR;
    cbtScalar dZ = dir.z() / pR;
    cbtScalar tpar = std::atan((pY * R_vert) / (pH * R_hor));
    cbtScalar sY = R_vert * std::sin(tpar);
    cbtScalar sH = R_hor * std::cos(tpar);
    cbtScalar sR = sH + R_offset;
    cbtScalar sX = dX * sR;
    cbtScalar sZ = dZ * sR;
    supVec.setValue(sX, sY, sZ);
    supVec -= dir * s;

    // support point on the top disc?
    if ((std::abs(Y_high) < R_vert) && (supVec.y() > sY_high)) {
        cbtScalar R_high = getCoreRadius(sY_high, s);
        supVecD.setValue(dX * R_high, sY_high, dZ * R_high);
        supVec = supVecD;
    }
    // support point on the bottom disc?
    if ((std::abs(Y_low) < R_vert) && (supVec.y() < sY_low)) {
        cbtScalar R_low = getCoreRadius(sY_low, s);
        supVecD.setValue(dX * R_low, sY_low, dZ * R_low);
        supVec = supVecD;
    }

    return supVec;
}

void cbtBarrelShape::batchedUnitVectorGetSupportingVertexWithoutMargin(const cbtVector3* vectors, cbtVector3* supportVerticesOut, int numVectors) const {
    for (int i = 0; i < numVectors; i++)
        supportVerticesOut[i] = localGetSupportingVertexWithoutMargin(vectors[i]);
}

void cbtBarrelShape::calculateLocalInertia(cbtScalar mass, cbtVector3& inertia) const {
    //// TODO
    // as an approximation, take the inertia of the box that bounds the barrel

    cbtTransform ident;
    ident.setIdentity();

    cbtVector3 halfExtents;

    halfExtents.setValue((R_hor + R_offset), (R_vert), (R_hor + R_offset));

    cbtScalar margin = CONVEX_DISTANCE_MARGIN;

    cbtScalar lx = cbtScalar(2.) * (halfExtents[0] + margin);
    cbtScalar ly = cbtScalar(2.) * (halfExtents[1] + margin);
    cbtScalar lz = cbtScalar(2.) * (halfExtents[2] + margin);
    const cbtScalar x2 = lx * lx;
    const cbtScalar y2 = ly * ly;
    const cbtScalar z2 = lz * lz;
    const cbtScalar scaledmass = mass * cbtScalar(.08333333);

    inertia[0] = scaledmass * (y2 + z2);
    inertia[1] = scaledmass * (x2 + z2);
    inertia[2] = scaledmass * (x2 + y2);
}

void cbtBarrelShape::getAabb(const cbtTransform& t, cbtVector3& aabbMin, cbtVector3& aabbMax) const {
    cbtVector3 halfExtents;
    halfExtents.setValue((R_hor + R_offset), (R_vert), (R_hor + R_offset));

    cbtMatrix3x3 abs_b = t.getBasis().absolute();
    cbtVector3 center = t.getOrigin();
    cbtVector3 extent = cbtVector3(abs_b[0].dot(halfExtents), abs_b[1].dot(halfExtents), abs_b[2].dot(halfExtents));
    extent += cbtVector3(m_offset, m_offset, m_offset);

    aabbMin = center - extent;
    aabbMax = center + extent;
}
