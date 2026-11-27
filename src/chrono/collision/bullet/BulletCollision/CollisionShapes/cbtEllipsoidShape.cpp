/*
***CHRONO***
*/

#include "cbtEllipsoidShape.h"

cbtEllipsoidShape::cbtEllipsoidShape(const cbtVector3& semiAxes, cbtScalar offset) : m_semiAxes(semiAxes), m_offset(offset) {
    m_shapeType = CUSTOM_CONVEX_SHAPE_TYPE;
}

cbtVector3 cbtEllipsoidShape::localGetSupportingVertexWithoutMargin(const cbtVector3& vec0) const {
    cbtVector3 vec = vec0;
    cbtScalar lenSqr = vec.length2();
    if (lenSqr < SIMD_EPSILON * SIMD_EPSILON)
        vec.setValue(1, 0, 0);
    else
        vec *= cbtScalar(1) / cbtSqrt(lenSqr);

    // support point of the ellipsoid x^T S^-2 x = 1 in the unit direction n: S*S*n / |S*n|
    cbtVector3 semi = m_semiAxes * m_localScaling;
    cbtVector3 sn = semi * vec;
    cbtScalar len = sn.length();
    cbtVector3 supVec = len > SIMD_EPSILON ? semi * sn / len : cbtVector3(0, 0, 0);

    // offset of the ellipsoid surface, minus the margin (which Bullet adds back along the same direction)
    return supVec + vec * (m_offset - getMargin());
}

void cbtEllipsoidShape::batchedUnitVectorGetSupportingVertexWithoutMargin(const cbtVector3* vectors, cbtVector3* supportVerticesOut, int numVectors) const {
    for (int i = 0; i < numVectors; i++)
        supportVerticesOut[i] = localGetSupportingVertexWithoutMargin(vectors[i]);
}

void cbtEllipsoidShape::calculateLocalInertia(cbtScalar mass, cbtVector3& inertia) const {
    // solid ellipsoid with the semi-axes inflated by the offset
    cbtVector3 semi = m_semiAxes * m_localScaling + cbtVector3(m_offset, m_offset, m_offset);
    cbtScalar a2 = semi.x() * semi.x();
    cbtScalar b2 = semi.y() * semi.y();
    cbtScalar c2 = semi.z() * semi.z();
    inertia.setValue(mass * (b2 + c2) / cbtScalar(5), mass * (a2 + c2) / cbtScalar(5), mass * (a2 + b2) / cbtScalar(5));
}

void cbtEllipsoidShape::getAabb(const cbtTransform& t, cbtVector3& aabbMin, cbtVector3& aabbMax) const {
    // extent of the rotated ellipsoid along a world axis e: |S * R^T * e|, plus the offset
    cbtVector3 semi = m_semiAxes * m_localScaling;
    const cbtMatrix3x3& basis = t.getBasis();
    cbtVector3 extent;
    for (int i = 0; i < 3; i++)
        extent[i] = (basis[i] * semi).length() + m_offset;

    cbtVector3 center = t.getOrigin();
    aabbMin = center - extent;
    aabbMax = center + extent;
}
