/*
***CHRONO***
*/

#ifndef BT_BARREL_SHAPE_H
#define BT_BARREL_SHAPE_H

#include "cbtConvexInternalShape.h"
#include "BulletCollision/BroadphaseCollision/cbtBroadphaseProxy.h"  // for the types

/// cbtBarrelShape represents a barrel-like shape (with vertical axis)
/// made by lathing an arc of an ellipse around the vertical Y axis.
/// The center of the ellipse is on Y=0 level, and it is offset by R_offset from
/// the Y axis in radial direction. The two radii of the ellipse are R_vert (for the
/// vertical direction, i.e. the axis parallel to Y) and R_hor (for the axis that
/// is perpendicular to Y). Also, the solid is clamped with two discs on the top and
/// the bottom, at levels Y_low and Y_high.
/// The shape can be inflated by a uniform outward offset of its surface. The collision margin is taken from the inside
/// of the shape: the shape including its margin is the (inflated) barrel.
class cbtBarrelShape : public cbtConvexInternalShape {
  private:
    cbtScalar Y_low;
    cbtScalar Y_high;
    cbtScalar R_vert;
    cbtScalar R_hor;
    cbtScalar R_offset;
    cbtScalar m_offset;

  public:
    cbtBarrelShape(cbtScalar sY_low, cbtScalar sY_high, cbtScalar sR_vert, cbtScalar sR_hor, cbtScalar sR_offset, cbtScalar offset = cbtScalar(0));

    /// CollisionShape Interface
    virtual void calculateLocalInertia(cbtScalar mass, cbtVector3& inertia) const;

    /// cbtConvexShape Interface
    virtual cbtVector3 localGetSupportingVertexWithoutMargin(const cbtVector3& vec) const;

    virtual void batchedUnitVectorGetSupportingVertexWithoutMargin(const cbtVector3* vectors, cbtVector3* supportVerticesOut, int numVectors) const;

    // virtual int	getShapeType() const { return BARREL_SHAPE_PROXYTYPE; }

    virtual const char* getName() const { return "BarrelShape"; }

    virtual void getAabb(const cbtTransform& t, cbtVector3& aabbMin, cbtVector3& aabbMax) const;

    cbtScalar getY_low() const { return Y_low; }
    cbtScalar getY_high() const { return Y_high; }
    cbtScalar getR_vert() const { return R_vert; }
    cbtScalar getR_hor() const { return R_hor; }
    cbtScalar getR_offset() const { return R_offset; }
    cbtScalar getOffset() const { return m_offset; }

  private:
    /// Radius of the lathed ellipse, offset inward by s, at the level y.
    cbtScalar getCoreRadius(cbtScalar y, cbtScalar s) const;
};

#endif
