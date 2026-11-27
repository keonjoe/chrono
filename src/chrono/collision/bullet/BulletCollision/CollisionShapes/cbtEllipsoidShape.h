/*
***CHRONO***
*/

#ifndef BT_ELLIPSOID_SHAPE_H
#define BT_ELLIPSOID_SHAPE_H

#include "cbtConvexInternalShape.h"
#include "BulletCollision/BroadphaseCollision/cbtBroadphaseProxy.h"  // for the types

/// cbtEllipsoidShape represents an ellipsoid centered at the origin, with semi-axes along the X, Y and Z axes,
/// optionally inflated by a uniform outward offset (i.e. the set of points within distance 'offset' from the
/// ellipsoid).
/// Unlike a unit cbtMultiSphereShape with non-uniform local scaling, whose support mapping returns S*v instead of
/// S*S*v/|S*v| (S = diag(semi-axes)) and therefore underestimates the extent of the ellipsoid in all directions that
/// are not aligned with an axis, the support mapping of this shape is exact.
/// The collision margin is taken from the inside of the shape: the shape including its margin is the (inflated)
/// ellipsoid. The core shape (without margin) is exact as long as the margin does not exceed the offset plus the
/// smallest radius of curvature of the ellipsoid surface, min(semi-axes)^2 / max(semi-axes).
ATTRIBUTE_ALIGNED16(class)
cbtEllipsoidShape : public cbtConvexInternalShape {
  private:
    cbtVector3 m_semiAxes;
    cbtScalar m_offset;

  public:
    BT_DECLARE_ALIGNED_ALLOCATOR();

    cbtEllipsoidShape(const cbtVector3& semiAxes, cbtScalar offset = cbtScalar(0));

    /// Semi-axes of the ellipsoid (without local scaling and offset).
    const cbtVector3& getSemiAxes() const { return m_semiAxes; }

    /// Uniform outward offset of the ellipsoid surface.
    cbtScalar getOffset() const { return m_offset; }

    /// CollisionShape Interface
    virtual void calculateLocalInertia(cbtScalar mass, cbtVector3 & inertia) const;

    /// cbtConvexShape Interface
    virtual cbtVector3 localGetSupportingVertexWithoutMargin(const cbtVector3& vec) const;

    virtual void batchedUnitVectorGetSupportingVertexWithoutMargin(const cbtVector3* vectors, cbtVector3* supportVerticesOut, int numVectors) const;

    virtual void getAabb(const cbtTransform& t, cbtVector3& aabbMin, cbtVector3& aabbMax) const;

    virtual const char* getName() const { return "Ellipsoid"; }
};

#endif
