/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  http://continuousphysics.com/Bullet/

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

///Specialized capsule-capsule collision algorithm has been added for Bullet 2.75 release to increase ragdoll performance
///If you experience problems with capsule-capsule collision, try to define BT_DISABLE_CAPSULE_CAPSULE_COLLIDER and report it in the Bullet forums
///with reproduction case
//#define BT_DISABLE_CAPSULE_CAPSULE_COLLIDER 1
//#define ZERO_MARGIN

#include "cbtConvexConvexAlgorithm.h"

//#include <stdio.h>
#include "BulletCollision/NarrowPhaseCollision/cbtDiscreteCollisionDetectorInterface.h"
#include "BulletCollision/BroadphaseCollision/cbtBroadphaseInterface.h"
#include "BulletCollision/CollisionDispatch/cbtCollisionObject.h"
#include "BulletCollision/CollisionShapes/cbtConvexShape.h"
#include "BulletCollision/CollisionShapes/cbtCapsuleShape.h"
#include "BulletCollision/CollisionShapes/cbtTriangleShape.h"
#include "BulletCollision/CollisionShapes/cbtConvexPolyhedron.h"

#include "BulletCollision/NarrowPhaseCollision/cbtGjkPairDetector.h"
#include "BulletCollision/BroadphaseCollision/cbtBroadphaseProxy.h"
#include "BulletCollision/CollisionDispatch/cbtCollisionDispatcher.h"
#include "BulletCollision/CollisionShapes/cbtBoxShape.h"
#include "BulletCollision/CollisionDispatch/cbtManifoldResult.h"

#include "BulletCollision/NarrowPhaseCollision/cbtConvexPenetrationDepthSolver.h"
#include "BulletCollision/NarrowPhaseCollision/cbtContinuousConvexCollision.h"
#include "BulletCollision/NarrowPhaseCollision/cbtSubSimplexConvexCast.h"
#include "BulletCollision/NarrowPhaseCollision/cbtGjkConvexCast.h"

#include "BulletCollision/NarrowPhaseCollision/cbtVoronoiSimplexSolver.h"
#include "BulletCollision/CollisionShapes/cbtSphereShape.h"
#include "BulletCollision/CollisionShapes/cbtEllipsoidShape.h" /* ***CHRONO*** */
#include "BulletCollision/CollisionShapes/cbtCylinderShape.h" /* ***CHRONO*** */
#include "BulletCollision/CollisionShapes/cbtRoundedBoxShape.h" /* ***CHRONO*** */
#include "BulletCollision/CollisionShapes/cbtRoundedCylinderShape.h" /* ***CHRONO*** */

#include "BulletCollision/NarrowPhaseCollision/cbtMinkowskiPenetrationDepthSolver.h"

#include "BulletCollision/NarrowPhaseCollision/cbtGjkEpa2.h"
#include "BulletCollision/NarrowPhaseCollision/cbtGjkEpaPenetrationDepthSolver.h"
#include "BulletCollision/NarrowPhaseCollision/cbtPolyhedralContactClipping.h"
#include "BulletCollision/CollisionDispatch/cbtCollisionObjectWrapper.h"

///////////

static SIMD_FORCE_INLINE void segmentsClosestPoints(
	cbtVector3& ptsVector,
	cbtVector3& offsetA,
	cbtVector3& offsetB,
	cbtScalar& tA, cbtScalar& tB,
	const cbtVector3& translation,
	const cbtVector3& dirA, cbtScalar hlenA,
	const cbtVector3& dirB, cbtScalar hlenB)
{
	// compute the parameters of the closest points on each line segment

	cbtScalar dirA_dot_dirB = cbtDot(dirA, dirB);
	cbtScalar dirA_dot_trans = cbtDot(dirA, translation);
	cbtScalar dirB_dot_trans = cbtDot(dirB, translation);

	cbtScalar denom = 1.0f - dirA_dot_dirB * dirA_dot_dirB;

	if (denom == 0.0f)
	{
		tA = 0.0f;
	}
	else
	{
		tA = (dirA_dot_trans - dirB_dot_trans * dirA_dot_dirB) / denom;
		if (tA < -hlenA)
			tA = -hlenA;
		else if (tA > hlenA)
			tA = hlenA;
	}

	tB = tA * dirA_dot_dirB - dirB_dot_trans;

	if (tB < -hlenB)
	{
		tB = -hlenB;
		tA = tB * dirA_dot_dirB + dirA_dot_trans;

		if (tA < -hlenA)
			tA = -hlenA;
		else if (tA > hlenA)
			tA = hlenA;
	}
	else if (tB > hlenB)
	{
		tB = hlenB;
		tA = tB * dirA_dot_dirB + dirA_dot_trans;

		if (tA < -hlenA)
			tA = -hlenA;
		else if (tA > hlenA)
			tA = hlenA;
	}

	// compute the closest points relative to segment centers.

	offsetA = dirA * tA;
	offsetB = dirB * tB;

	ptsVector = translation - offsetA + offsetB;
}

static SIMD_FORCE_INLINE cbtScalar capsuleCapsuleDistance(
	cbtVector3& normalOnB,
	cbtVector3& pointOnB,
	cbtScalar capsuleLengthA,
	cbtScalar capsuleRadiusA,
	cbtScalar capsuleLengthB,
	cbtScalar capsuleRadiusB,
	int capsuleAxisA,
	int capsuleAxisB,
	const cbtTransform& transformA,
	const cbtTransform& transformB,
	cbtScalar distanceThreshold)
{
	cbtVector3 directionA = transformA.getBasis().getColumn(capsuleAxisA);
	cbtVector3 translationA = transformA.getOrigin();
	cbtVector3 directionB = transformB.getBasis().getColumn(capsuleAxisB);
	cbtVector3 translationB = transformB.getOrigin();

	// translation between centers

	cbtVector3 translation = translationB - translationA;

	// compute the closest points of the capsule line segments

	cbtVector3 ptsVector;  // the vector between the closest points

	cbtVector3 offsetA, offsetB;  // offsets from segment centers to their closest points
	cbtScalar tA, tB;             // parameters on line segment

	segmentsClosestPoints(ptsVector, offsetA, offsetB, tA, tB, translation,
						  directionA, capsuleLengthA, directionB, capsuleLengthB);

	cbtScalar distance = ptsVector.length() - capsuleRadiusA - capsuleRadiusB;

	if (distance > distanceThreshold)
		return distance;

	cbtScalar lenSqr = ptsVector.length2();
	if (lenSqr <= (SIMD_EPSILON * SIMD_EPSILON))
	{
		//degenerate case where 2 capsules are likely at the same location: take a vector tangential to 'directionA'
		cbtVector3 q;
		cbtPlaneSpace1(directionA, normalOnB, q);
	}
	else
	{
		// compute the contact normal
		normalOnB = ptsVector * -cbtRecipSqrt(lenSqr);
	}
	pointOnB = transformB.getOrigin() + offsetB + normalOnB * capsuleRadiusB;

	return distance;
}

//////////

cbtConvexConvexAlgorithm::CreateFunc::CreateFunc(cbtConvexPenetrationDepthSolver* pdSolver)
{
	/* ***CHRONO*** GJK returns one point per query. Without perturbation a flat face resting on a flat face (convex
	hulls, cylinder and cone bases, boxes on non-box shapes) is supported by a single point that wanders over the face;
	the perturbed queries add the corners of the contact region (as cbtDefaultCollisionConfiguration::
	setConvexConvexMultipointIterations does with its default arguments). */
	m_numPerturbationIterations = 3;
	m_minimumPointsPerturbationThreshold = 3;
	m_pdSolver = pdSolver;
}

cbtConvexConvexAlgorithm::CreateFunc::~CreateFunc()
{
}

cbtConvexConvexAlgorithm::cbtConvexConvexAlgorithm(cbtPersistentManifold* mf, const cbtCollisionAlgorithmConstructionInfo& ci, const cbtCollisionObjectWrapper* body0Wrap, const cbtCollisionObjectWrapper* body1Wrap, cbtConvexPenetrationDepthSolver* pdSolver, int numPerturbationIterations, int minimumPointsPerturbationThreshold)
	: cbtActivatingCollisionAlgorithm(ci, body0Wrap, body1Wrap),
	  m_pdSolver(pdSolver),
	  m_ownManifold(false),
	  m_manifoldPtr(mf),
	  m_lowLevelOfDetail(false),
#ifdef USE_SEPDISTANCE_UTIL2
	  m_sepDistance((static_cast<cbtConvexShape*>(body0->getCollisionShape()))->getAngularMotionDisc(),
					(static_cast<cbtConvexShape*>(body1->getCollisionShape()))->getAngularMotionDisc()),
#endif
	  m_numPerturbationIterations(numPerturbationIterations),
	  m_minimumPointsPerturbationThreshold(minimumPointsPerturbationThreshold)
{
	(void)body0Wrap;
	(void)body1Wrap;
}

cbtConvexConvexAlgorithm::~cbtConvexConvexAlgorithm()
{
	if (m_ownManifold)
	{
		if (m_manifoldPtr)
			m_dispatcher->releaseManifold(m_manifoldPtr);
	}
}

void cbtConvexConvexAlgorithm ::setLowLevelOfDetail(bool useLowLevel)
{
	m_lowLevelOfDetail = useLowLevel;
}

struct cbtPerturbedContactResult : public cbtManifoldResult
{
	cbtManifoldResult* m_originalManifoldResult;
	cbtTransform m_transformA;
	cbtTransform m_transformB;
	cbtTransform m_unPerturbedTransform;
	bool m_perturbA;
	cbtIDebugDraw* m_debugDrawer;
	/* ***CHRONO*** collision margin of the perturbed shape; normal of the unperturbed contact and the tolerance on the
	cosine of the angle between it and the normal of a perturbed contact */
	cbtScalar m_margin;
	cbtVector3 m_normal;
	cbtScalar m_normalTolerance;

	cbtPerturbedContactResult(cbtManifoldResult* originalResult, const cbtTransform& transformA, const cbtTransform& transformB, const cbtTransform& unPerturbedTransform, bool perturbA, cbtIDebugDraw* debugDrawer,
							  cbtScalar margin = cbtScalar(0), const cbtVector3& normal = cbtVector3(0, 0, 0), cbtScalar normalTolerance = cbtScalar(2))
		: m_originalManifoldResult(originalResult),
		  m_transformA(transformA),
		  m_transformB(transformB),
		  m_unPerturbedTransform(unPerturbedTransform),
		  m_perturbA(perturbA),
		  m_debugDrawer(debugDrawer),
		  m_margin(margin),
		  m_normal(normal),
		  m_normalTolerance(normalTolerance)
	{
	}
	virtual ~cbtPerturbedContactResult()
	{
	}

	virtual void addContactPoint(const cbtVector3& normalOnBInWorld, const cbtVector3& pointInWorld, cbtScalar orgDepth)
	{
		/* ***CHRONO*** a perturbed contact is a valid extra point only where a feature of the perturbed shape touches a flat
		face of the other one, i.e. with the normal of the unperturbed contact (corners of a face resting on a face, ends
		of an edge or a cylinder resting on a face). On curved surfaces the normal turns with the perturbation and the
		point is not a contact point of the unperturbed configuration. */
		if (cbtFabs(normalOnBInWorld.dot(m_normal)) < cbtScalar(1) - m_normalTolerance)
			return;

		cbtVector3 endPt, startPt;
		cbtScalar newDepth;
		cbtVector3 newNormal;

		/* ***CHRONO*** the witness point of the perturbed shape lies on its margin, offset from the core point along the
		normal of the perturbed configuration: move the core point back to the unperturbed configuration and apply the
		margin along the normal there (otherwise the depth of a point on a flat face is off by margin * (1 - cos(angle))) */
		if (m_perturbA)
		{
			cbtVector3 endPtOrg = pointInWorld + normalOnBInWorld * orgDepth;
			endPt = (m_unPerturbedTransform * m_transformA.inverse())(endPtOrg + normalOnBInWorld * m_margin) - normalOnBInWorld * m_margin;
			newDepth = (endPt - pointInWorld).dot(normalOnBInWorld);
			startPt = endPt - normalOnBInWorld * newDepth;
		}
		else
		{
			endPt = pointInWorld + normalOnBInWorld * orgDepth;
			startPt = (m_unPerturbedTransform * m_transformB.inverse())(pointInWorld - normalOnBInWorld * m_margin) + normalOnBInWorld * m_margin;
			newDepth = (endPt - startPt).dot(normalOnBInWorld);
		}

//#define DEBUG_CONTACTS 1
#ifdef DEBUG_CONTACTS
		m_debugDrawer->drawLine(startPt, endPt, cbtVector3(1, 0, 0));
		m_debugDrawer->drawSphere(startPt, 0.05, cbtVector3(0, 1, 0));
		m_debugDrawer->drawSphere(endPt, 0.05, cbtVector3(0, 0, 1));
#endif  //DEBUG_CONTACTS

		/* ***CHRONO*** a perturbed point close to a point already in the manifold (e.g. the unperturbed point of a curved
		contact) would replace it: keep only the points that add a new location (the corners of a flat contact region) */
		const cbtPersistentManifold* manifold = m_originalManifoldResult->getPersistentManifold();
		cbtScalar threshold2 = manifold->getContactBreakingThreshold() * manifold->getContactBreakingThreshold();
		cbtVector3 pointOnA = startPt + normalOnBInWorld * newDepth;
		for (int i = 0; i < manifold->getNumContacts(); i++)
		{
			const cbtManifoldPoint& pt = manifold->getContactPoint(i);
			if ((pt.getPositionWorldOnB() - startPt).length2() < threshold2 || (pt.getPositionWorldOnA() - pointOnA).length2() < threshold2)
				return;
		}

		m_originalManifoldResult->addContactPoint(normalOnBInWorld, startPt, newDepth);
	}
};

extern cbtScalar gContactBreakingThreshold;

/* ***CHRONO*** */
static bool isStrictlyConvex(const cbtConvexShape* shape)
{
	return shape->getShapeType() == SPHERE_SHAPE_PROXYTYPE || dynamic_cast<const cbtEllipsoidShape*>(shape) != nullptr;
}

/* ***CHRONO*** triangles of a mesh (rigid or FE contact surface), segments of a beam and points of a node cloud are
elements of a discretised surface, line or cloud: the corners that perturbing such an element adds are mesh nodes,
which the neighbouring elements already support. The manifold of a triangle is cleared at each step (no contact
persistence), so perturbing the triangles would add up to three points per triangle in every step, multiplying the
contact stiffness of the surface. The other shape of the pair (e.g. a box resting on a large triangle) is still
perturbed. */
static bool isMeshElement(const cbtConvexShape* shape)
{
	int type = shape->getShapeType();
	return type == CH_TRIANGLE_SHAPE_PROXYTYPE || type == SEGMENT_SHAPE_PROXYTYPE || type == POINT_SHAPE_PROXYTYPE;
}

/* ***CHRONO*** Exact contact of a sphere with a box, a cylinder, a rounded box or a rounded cylinder. These shapes are a
core box or cylinder (the half extents without margin, less the sweep radius of the rounded shapes) swept by a sphere
(the margin, plus the sweep radius), so the closest point to the sphere center follows from clamping the center to the
core. GJK reaches the same point only through simplices with support points on the far rims or corners of the core,
whose closest point loses accuracy in proportion to the length of the shape (e.g. 1e-6 m of the contact position on a
10 m cylinder in double precision), which produces spurious torques on a sphere rolling along the shape. */
/* ***CHRONO*** Contact normal on the lateral side of a cylinder, a rounded cylinder or a capsule. GJK finds the normal
as the direction of the separation vector of the core shapes, whose length is only the sum of the margins; on a long
shape its error is a sizable fraction of that length (the simplex has vertices at both far ends), which tilts the
normal, also along the axis, and produces spurious forces along the shape. On the lateral side the normal is radial,
so it follows from the witness point alone, with an error of (witness point error) / radius instead of
(witness point error) / (sum of the margins), and without axial component. Not used against polyhedra and mesh
elements, whose flat faces give GJK an exact normal. */
static bool lateralCore(const cbtConvexShape* shape, int& axis, cbtScalar& halfLength)
{
	switch (shape->getShapeType())
	{
		case CYLINDER_SHAPE_PROXYTYPE:
			axis = static_cast<const cbtCylinderShape*>(shape)->getUpAxis();
			halfLength = static_cast<const cbtCylinderShape*>(shape)->getHalfExtentsWithoutMargin()[axis];
			return true;
		case ROUNDEDCYL_SHAPE_PROXYTYPE:
		{
			const cbtRoundedCylinderShape* rcyl = static_cast<const cbtRoundedCylinderShape*>(shape);
			axis = 2;
			halfLength = rcyl->getHalfExtentsWithoutMargin()[2] - rcyl->getSphereRadius();
			return true;
		}
		case CAPSULE_SHAPE_PROXYTYPE:
			axis = static_cast<const cbtCapsuleShape*>(shape)->getUpAxis();
			halfLength = static_cast<const cbtCapsuleShape*>(shape)->getHalfHeight();
			return true;
		default:
			return false;
	}
}

// Outward radial direction and axis (world) at the world point p on the lateral side of the shape; false if p is
// beyond the lateral side (cap, rim) or on the axis.
static bool lateralNormal(const cbtConvexShape* shape, const cbtTransform& trans, const cbtVector3& p, cbtVector3& normal, cbtVector3& axisDir)
{
	int axis;
	cbtScalar halfLength;
	if (!lateralCore(shape, axis, halfLength))
		return false;
	cbtVector3 r = trans.invXform(p);
	if (cbtFabs(r[axis]) >= halfLength)
		return false;
	r[axis] = cbtScalar(0);
	cbtScalar radius = r.length();
	if (radius <= shape->getMargin())
		return false;
	normal = trans.getBasis() * (r / radius);
	axisDir = trans.getBasis().getColumn(axis);
	return true;
}

// The normal without its component along the unit vector 'axis'
static cbtVector3 removeAxialComponent(const cbtVector3& normal, const cbtVector3& axis)
{
	cbtVector3 n = normal - axis * normal.dot(axis);
	cbtScalar len = n.length();
	return len > cbtScalar(0.5) ? n / len : normal;
}

struct cbtLateralNormalResult : public cbtDiscreteCollisionDetectorInterface::Result
{
	cbtDiscreteCollisionDetectorInterface::Result* m_original;
	const cbtConvexShape* m_shapeA;
	const cbtConvexShape* m_shapeB;
	cbtTransform m_transA;
	cbtTransform m_transB;
	bool m_corrected;
	cbtVector3 m_normal;

	cbtLateralNormalResult(cbtDiscreteCollisionDetectorInterface::Result* original, const cbtConvexShape* shapeA, const cbtConvexShape* shapeB,
						   const cbtTransform& transA, const cbtTransform& transB)
		: m_original(original), m_shapeA(shapeA), m_shapeB(shapeB), m_transA(transA), m_transB(transB), m_corrected(false)
	{
	}

	static bool applies(const cbtConvexShape* shapeA, const cbtConvexShape* shapeB)
	{
		int axis;
		cbtScalar halfLength;
		bool lateralA = lateralCore(shapeA, axis, halfLength) && !shapeB->isPolyhedral() && !isMeshElement(shapeB);
		bool lateralB = lateralCore(shapeB, axis, halfLength) && !shapeA->isPolyhedral() && !isMeshElement(shapeA);
		return lateralA || lateralB;
	}

	virtual void setShapeIdentifiersA(int partId0, int index0) { m_original->setShapeIdentifiersA(partId0, index0); }
	virtual void setShapeIdentifiersB(int partId1, int index1) { m_original->setShapeIdentifiersB(partId1, index1); }
	virtual void addContactPoint(const cbtVector3& normalOnBInWorld, const cbtVector3& pointInWorld, cbtScalar depth)
	{
		// normalOnBInWorld points from B to A; the witness point on A is pointInWorld + normalOnBInWorld * depth.
		// On the lateral side of a shape the normal has no axial component; if both witness points are on lateral
		// sides of shapes with crossed axes, the normal is along the cross product of the axes. Against a strictly
		// convex shape the radial direction through the witness point is used.
		cbtVector3 normal = normalOnBInWorld;
		cbtVector3 nA, nB, axisA, axisB;
		bool lateralA = !m_shapeB->isPolyhedral() && !isMeshElement(m_shapeB) && lateralNormal(m_shapeA, m_transA, pointInWorld + normalOnBInWorld * depth, nA, axisA);
		bool lateralB = !m_shapeA->isPolyhedral() && !isMeshElement(m_shapeA) && lateralNormal(m_shapeB, m_transB, pointInWorld, nB, axisB);
		if (lateralA && lateralB)
		{
			cbtVector3 c = axisA.cross(axisB);
			cbtScalar len = c.length();
			if (len > cbtScalar(1e-3))
			{
				c /= len;
				normal = c.dot(normalOnBInWorld) < 0 ? -c : c;
			}
			else
				normal = removeAxialComponent(normalOnBInWorld, axisA);
		}
		else if (lateralA)
			normal = (isStrictlyConvex(m_shapeB) && -nA.dot(normalOnBInWorld) > 0) ? -nA : removeAxialComponent(normalOnBInWorld, axisA);
		else if (lateralB)
			normal = (isStrictlyConvex(m_shapeA) && nB.dot(normalOnBInWorld) > 0) ? nB : removeAxialComponent(normalOnBInWorld, axisB);
		if (lateralA || lateralB)
		{
			m_corrected = true;
			m_normal = normal;
		}
		m_original->addContactPoint(normal, pointInWorld, depth);
	}
};

static bool sphereCorePrimitive(const cbtConvexShape* shape, int& kind, cbtVector3& halfExtents, int& axis, cbtScalar& rounding)
{
	switch (shape->getShapeType())
	{
		case BOX_SHAPE_PROXYTYPE:
			kind = 0;
			halfExtents = static_cast<const cbtBoxShape*>(shape)->getHalfExtentsWithoutMargin();
			rounding = shape->getMargin();
			return true;
		case CYLINDER_SHAPE_PROXYTYPE:
			kind = 1;
			halfExtents = static_cast<const cbtCylinderShape*>(shape)->getHalfExtentsWithoutMargin();
			axis = static_cast<const cbtCylinderShape*>(shape)->getUpAxis();
			rounding = shape->getMargin();
			return true;
		case ROUNDEDBOX_SHAPE_PROXYTYPE:
		{
			const cbtRoundedBoxShape* rbox = static_cast<const cbtRoundedBoxShape*>(shape);
			cbtScalar sr = rbox->getSphereRadius();
			kind = 0;
			halfExtents = rbox->getHalfExtentsWithoutMargin() - cbtVector3(sr, sr, sr);
			rounding = shape->getMargin() + sr;
			return true;
		}
		case ROUNDEDCYL_SHAPE_PROXYTYPE:
		{
			const cbtRoundedCylinderShape* rcyl = static_cast<const cbtRoundedCylinderShape*>(shape);
			cbtScalar sr = rcyl->getSphereRadius();
			kind = 1;
			halfExtents = rcyl->getHalfExtentsWithoutMargin() - cbtVector3(sr, sr, sr);
			axis = 2;
			rounding = shape->getMargin() + sr;
			return true;
		}
		default:
			return false;
	}
}

// Signed distance from the point p to the surface of a core box (kind 0) or cylinder (kind 1, axis 'axis'), in the
// frame of the core; q is the closest point on the core surface and n the outward unit normal there.
static cbtScalar closestPointOnCore(int kind, const cbtVector3& h, int axis, const cbtVector3& p, cbtVector3& q, cbtVector3& n)
{
	if (kind == 0)
	{
		q.setValue(cbtClamped(p.x(), -h.x(), h.x()), cbtClamped(p.y(), -h.y(), h.y()), cbtClamped(p.z(), -h.z(), h.z()));
		cbtVector3 d = p - q;
		cbtScalar len = d.length();
		if (len > SIMD_EPSILON * h.length())
		{
			n = d / len;
			return len;
		}
		// inside: nearest face
		int i = 0;
		cbtScalar best = h[0] - cbtFabs(p[0]);
		for (int j = 1; j < 3; j++)
		{
			if (h[j] - cbtFabs(p[j]) < best)
			{
				best = h[j] - cbtFabs(p[j]);
				i = j;
			}
		}
		n.setValue(0, 0, 0);
		n[i] = p[i] < 0 ? cbtScalar(-1) : cbtScalar(1);
		q = p;
		q[i] = n[i] * h[i];
		return -best;
	}

	int i1 = (axis + 1) % 3, i2 = (axis + 2) % 3;
	cbtScalar r = h[i1];
	cbtScalar hh = h[axis];
	cbtScalar rl = cbtSqrt(p[i1] * p[i1] + p[i2] * p[i2]);
	if (rl > r || cbtFabs(p[axis]) > hh)
	{
		q = p;
		if (rl > r)
		{
			q[i1] = p[i1] * r / rl;
			q[i2] = p[i2] * r / rl;
		}
		q[axis] = cbtClamped(p[axis], -hh, hh);
		cbtVector3 d = p - q;
		cbtScalar len = d.length();
		if (len > SIMD_EPSILON * (r + hh))
		{
			n = d / len;
			return len;
		}
	}
	// inside (or on the surface): nearest of the lateral surface and the end faces
	cbtScalar dSide = r - rl, dCap = hh - cbtFabs(p[axis]);
	n.setValue(0, 0, 0);
	q = p;
	if (dSide < dCap)
	{
		if (rl > SIMD_EPSILON * r)
		{
			n[i1] = p[i1] / rl;
			n[i2] = p[i2] / rl;
		}
		else
			n[i1] = 1;
		q[i1] = n[i1] * r;
		q[i2] = n[i2] * r;
		return -dSide;
	}
	n[axis] = p[axis] < 0 ? cbtScalar(-1) : cbtScalar(1);
	q[axis] = n[axis] * hh;
	return -dCap;
}

//
// Convex-Convex collision algorithm
//
void cbtConvexConvexAlgorithm ::processCollision(const cbtCollisionObjectWrapper* body0Wrap, const cbtCollisionObjectWrapper* body1Wrap, const cbtDispatcherInfo& dispatchInfo, cbtManifoldResult* resultOut)
{
	if (!m_manifoldPtr)
	{
		//swapped?
		m_manifoldPtr = m_dispatcher->getNewManifold(body0Wrap->getCollisionObject(), body1Wrap->getCollisionObject());
		m_ownManifold = true;
	}
	resultOut->setPersistentManifold(m_manifoldPtr);

	//comment-out next line to test multi-contact generation
	//resultOut->getPersistentManifold()->clearManifold();

	const cbtConvexShape* min0 = static_cast<const cbtConvexShape*>(body0Wrap->getCollisionShape());
	const cbtConvexShape* min1 = static_cast<const cbtConvexShape*>(body1Wrap->getCollisionShape());

	/* ***CHRONO*** disable contact persistence for deformable triangle meshes vs convexes, 
	because these do not simply move/rotate with affine transformation as other primitives. */
	if((min0->getShapeType() == CH_TRIANGLE_SHAPE_PROXYTYPE) ||  (min1->getShapeType() == CH_TRIANGLE_SHAPE_PROXYTYPE)) 
		resultOut->getPersistentManifold()->clearManifold();

	/* ***CHRONO*** a strictly convex shape (sphere, ellipsoid) touches a convex shape in a single point. Points kept
	in the manifold from earlier steps are refreshed with their old normal and, as the contact point moves on the curved
	surface, report penetrations that do not exist: keep only the current point. */
	const bool strictlyConvexPair = isStrictlyConvex(min0) || isStrictlyConvex(min1);
	if (m_ownManifold && strictlyConvexPair)
		resultOut->getPersistentManifold()->clearManifold();

	cbtVector3 normalOnB;
	cbtVector3 pointOnBWorld;

	/* ***CHRONO*** exact contact of a sphere with a box, a cylinder, a rounded box or a rounded cylinder */
	{
		bool sphere0 = min0->getShapeType() == SPHERE_SHAPE_PROXYTYPE;
		bool sphere1 = min1->getShapeType() == SPHERE_SHAPE_PROXYTYPE;
		int kind = 0, axis = 2;
		cbtVector3 halfExtents;
		cbtScalar rounding = 0;
		if ((sphere0 || sphere1) && sphereCorePrimitive(sphere0 ? min1 : min0, kind, halfExtents, axis, rounding))
		{
			const cbtTransform& sphereTrans = (sphere0 ? body0Wrap : body1Wrap)->getWorldTransform();
			const cbtTransform& coreTrans = (sphere0 ? body1Wrap : body0Wrap)->getWorldTransform();
			cbtScalar radius = static_cast<const cbtSphereShape*>(sphere0 ? min0 : min1)->getRadius();
			cbtVector3 q, n;
			cbtScalar coreDist = closestPointOnCore(kind, halfExtents, axis, coreTrans.invXform(sphereTrans.getOrigin()), q, n);
			cbtScalar dist = coreDist - rounding - radius;
			cbtVector3 normal = coreTrans.getBasis() * n;  // from the core shape towards the sphere center
			cbtScalar threshold = m_manifoldPtr->getContactBreakingThreshold() + resultOut->m_closestPointDistanceThreshold;
			if (dist < threshold)
			{
				if (sphere0)
					resultOut->addContactPoint(normal, coreTrans(q) + normal * rounding, dist);
				else
					resultOut->addContactPoint(-normal, sphereTrans.getOrigin() - normal * radius, dist);
			}
			if (m_ownManifold)
				resultOut->refreshContactPoints();
			return;
		}
	}

#ifndef BT_DISABLE_CAPSULE_CAPSULE_COLLIDER
	if ((min0->getShapeType() == CAPSULE_SHAPE_PROXYTYPE) && (min1->getShapeType() == CAPSULE_SHAPE_PROXYTYPE))
	{
		//m_manifoldPtr->clearManifold();

		cbtCapsuleShape* capsuleA = (cbtCapsuleShape*)min0;
		cbtCapsuleShape* capsuleB = (cbtCapsuleShape*)min1;

		cbtScalar threshold = m_manifoldPtr->getContactBreakingThreshold()+ resultOut->m_closestPointDistanceThreshold;

		cbtScalar dist = capsuleCapsuleDistance(normalOnB, pointOnBWorld, capsuleA->getHalfHeight(), capsuleA->getRadius(),
											   capsuleB->getHalfHeight(), capsuleB->getRadius(), capsuleA->getUpAxis(), capsuleB->getUpAxis(),
											   body0Wrap->getWorldTransform(), body1Wrap->getWorldTransform(), threshold);

		if (dist < threshold)
		{
			cbtAssert(normalOnB.length2() >= (SIMD_EPSILON * SIMD_EPSILON));
			resultOut->addContactPoint(normalOnB, pointOnBWorld, dist);
		}
		resultOut->refreshContactPoints();
		return;
	}

	if ((min0->getShapeType() == CAPSULE_SHAPE_PROXYTYPE) && (min1->getShapeType() == SPHERE_SHAPE_PROXYTYPE))
	{
		//m_manifoldPtr->clearManifold();

		cbtCapsuleShape* capsuleA = (cbtCapsuleShape*)min0;
		cbtSphereShape* capsuleB = (cbtSphereShape*)min1;

		cbtScalar threshold = m_manifoldPtr->getContactBreakingThreshold()+ resultOut->m_closestPointDistanceThreshold;

		cbtScalar dist = capsuleCapsuleDistance(normalOnB, pointOnBWorld, capsuleA->getHalfHeight(), capsuleA->getRadius(),
											   0., capsuleB->getRadius(), capsuleA->getUpAxis(), 1,
											   body0Wrap->getWorldTransform(), body1Wrap->getWorldTransform(), threshold);

		if (dist < threshold)
		{
			cbtAssert(normalOnB.length2() >= (SIMD_EPSILON * SIMD_EPSILON));
			resultOut->addContactPoint(normalOnB, pointOnBWorld, dist);
		}
		resultOut->refreshContactPoints();
		return;
	}

	if ((min0->getShapeType() == SPHERE_SHAPE_PROXYTYPE) && (min1->getShapeType() == CAPSULE_SHAPE_PROXYTYPE))
	{
		//m_manifoldPtr->clearManifold();

		cbtSphereShape* capsuleA = (cbtSphereShape*)min0;
		cbtCapsuleShape* capsuleB = (cbtCapsuleShape*)min1;

		cbtScalar threshold = m_manifoldPtr->getContactBreakingThreshold()+ resultOut->m_closestPointDistanceThreshold;

		cbtScalar dist = capsuleCapsuleDistance(normalOnB, pointOnBWorld, 0., capsuleA->getRadius(),
											   capsuleB->getHalfHeight(), capsuleB->getRadius(), 1, capsuleB->getUpAxis(),
											   body0Wrap->getWorldTransform(), body1Wrap->getWorldTransform(), threshold);

		if (dist < threshold)
		{
			cbtAssert(normalOnB.length2() >= (SIMD_EPSILON * SIMD_EPSILON));
			resultOut->addContactPoint(normalOnB, pointOnBWorld, dist);
		}
		resultOut->refreshContactPoints();
		return;
	}
#endif  //BT_DISABLE_CAPSULE_CAPSULE_COLLIDER

#ifdef USE_SEPDISTANCE_UTIL2
	if (dispatchInfo.m_useConvexConservativeDistanceUtil)
	{
		m_sepDistance.updateSeparatingDistance(body0->getWorldTransform(), body1->getWorldTransform());
	}

	if (!dispatchInfo.m_useConvexConservativeDistanceUtil || m_sepDistance.getConservativeSeparatingDistance() <= 0.f)
#endif  //USE_SEPDISTANCE_UTIL2

	{
		cbtGjkPairDetector::ClosestPointInput input;
		cbtVoronoiSimplexSolver simplexSolver;
		/* ***CHRONO*** */
		// The default threshold for detecting a repeated support point (a squared distance of 1e-4, i.e. 1 cm in
		// single precision) stops GJK early on small curved shapes (ellipsoids, cylinders, capsules, ...), with
		// distance errors of a sizable fraction of a millimeter. Scale it with the collision margins, which are
		// proportional to the size of the shapes. Pairs of polyhedra are left unchanged: their support points
		// are vertices, which repeat exactly.
		if (!min0->isPolyhedral() || !min1->isPolyhedral())
		{
			cbtScalar eqLen = cbtScalar(0.01) * cbtMin(min0->getMargin(), min1->getMargin());
			simplexSolver.setEqualVertexThreshold(cbtMin(simplexSolver.getEqualVertexThreshold(), eqLen * eqLen));
		}
		cbtGjkPairDetector gjkPairDetector(min0, min1, &simplexSolver, m_pdSolver);
		//TODO: if (dispatchInfo.m_useContinuous)
		gjkPairDetector.setMinkowskiA(min0);
		gjkPairDetector.setMinkowskiB(min1);

#ifdef USE_SEPDISTANCE_UTIL2
		if (dispatchInfo.m_useConvexConservativeDistanceUtil)
		{
			input.m_maximumDistanceSquared = BT_LARGE_FLOAT;
		}
		else
#endif  //USE_SEPDISTANCE_UTIL2
		{
			//if (dispatchInfo.m_convexMaxDistanceUseCPT)
			//{
			//	input.m_maximumDistanceSquared = min0->getMargin() + min1->getMargin() + m_manifoldPtr->getContactProcessingThreshold();
			//} else
			//{
			input.m_maximumDistanceSquared = min0->getMargin() + min1->getMargin() + m_manifoldPtr->getContactBreakingThreshold() + resultOut->m_closestPointDistanceThreshold;
			//		}

			input.m_maximumDistanceSquared *= input.m_maximumDistanceSquared;
		}

		input.m_transformA = body0Wrap->getWorldTransform();
		input.m_transformB = body1Wrap->getWorldTransform();

#ifdef USE_SEPDISTANCE_UTIL2
		cbtScalar sepDist = 0.f;
		if (dispatchInfo.m_useConvexConservativeDistanceUtil)
		{
			sepDist = gjkPairDetector.getCachedSeparatingDistance();
			if (sepDist > SIMD_EPSILON)
			{
				sepDist += dispatchInfo.m_convexConservativeDistanceThreshold;
				//now perturbe directions to get multiple contact points
			}
		}
#endif  //USE_SEPDISTANCE_UTIL2

		if (min0->isPolyhedral() && min1->isPolyhedral())
		{
			struct cbtDummyResult : public cbtDiscreteCollisionDetectorInterface::Result
			{
				cbtVector3 m_normalOnBInWorld;
				cbtVector3 m_pointInWorld;
				cbtScalar m_depth;
				bool m_hasContact;

				cbtDummyResult()
					: m_hasContact(false)
				{
				}

				virtual void setShapeIdentifiersA(int partId0, int index0) {}
				virtual void setShapeIdentifiersB(int partId1, int index1) {}
				virtual void addContactPoint(const cbtVector3& normalOnBInWorld, const cbtVector3& pointInWorld, cbtScalar depth)
				{
					m_hasContact = true;
					m_normalOnBInWorld = normalOnBInWorld;
					m_pointInWorld = pointInWorld;
					m_depth = depth;
				}
			};

			struct cbtWithoutMarginResult : public cbtDiscreteCollisionDetectorInterface::Result
			{
				cbtDiscreteCollisionDetectorInterface::Result* m_originalResult;
				cbtVector3 m_reportedNormalOnWorld;
				cbtScalar m_marginOnA;
				cbtScalar m_marginOnB;
				cbtScalar m_reportedDistance;

				bool m_foundResult;
				cbtWithoutMarginResult(cbtDiscreteCollisionDetectorInterface::Result* result, cbtScalar marginOnA, cbtScalar marginOnB)
					: m_originalResult(result),
					  m_marginOnA(marginOnA),
					  m_marginOnB(marginOnB),
					  m_foundResult(false)
				{
				}

				virtual void setShapeIdentifiersA(int partId0, int index0) {}
				virtual void setShapeIdentifiersB(int partId1, int index1) {}
				virtual void addContactPoint(const cbtVector3& normalOnBInWorld, const cbtVector3& pointInWorldOrg, cbtScalar depthOrg)
				{
					m_reportedDistance = depthOrg;
					m_reportedNormalOnWorld = normalOnBInWorld;

					cbtVector3 adjustedPointB = pointInWorldOrg - normalOnBInWorld * m_marginOnB;
					m_reportedDistance = depthOrg + (m_marginOnA + m_marginOnB);
					if (m_reportedDistance < 0.f)
					{
						m_foundResult = true;
					}
					m_originalResult->addContactPoint(normalOnBInWorld, adjustedPointB, m_reportedDistance);
				}
			};

			cbtDummyResult dummy;

			///cbtBoxShape is an exception: its vertices are created WITH margin so don't subtract it

			cbtScalar min0Margin = min0->getShapeType() == BOX_SHAPE_PROXYTYPE ? 0.f : min0->getMargin();
			cbtScalar min1Margin = min1->getShapeType() == BOX_SHAPE_PROXYTYPE ? 0.f : min1->getMargin();

			cbtWithoutMarginResult withoutMargin(resultOut, min0Margin, min1Margin);

			cbtPolyhedralConvexShape* polyhedronA = (cbtPolyhedralConvexShape*)min0;
			cbtPolyhedralConvexShape* polyhedronB = (cbtPolyhedralConvexShape*)min1;
			if (polyhedronA->getConvexPolyhedron() && polyhedronB->getConvexPolyhedron())
			{
				cbtScalar threshold = m_manifoldPtr->getContactBreakingThreshold()+ resultOut->m_closestPointDistanceThreshold;

				cbtScalar minDist = -1e30f;
				cbtVector3 sepNormalWorldSpace;
				bool foundSepAxis = true;

				if (dispatchInfo.m_enableSatConvex)
				{
					foundSepAxis = cbtPolyhedralContactClipping::findSeparatingAxis(
						*polyhedronA->getConvexPolyhedron(), *polyhedronB->getConvexPolyhedron(),
						body0Wrap->getWorldTransform(),
						body1Wrap->getWorldTransform(),
						sepNormalWorldSpace, *resultOut);
				}
				else
				{
#ifdef ZERO_MARGIN
					gjkPairDetector.setIgnoreMargin(true);
					gjkPairDetector.getClosestPoints(input, *resultOut, dispatchInfo.m_debugDraw);
#else

					gjkPairDetector.getClosestPoints(input, withoutMargin, dispatchInfo.m_debugDraw);
					//gjkPairDetector.getClosestPoints(input,dummy,dispatchInfo.m_debugDraw);
#endif  //ZERO_MARGIN                                                    \
	//cbtScalar l2 = gjkPairDetector.getCachedSeparatingAxis().length2(); \
	//if (l2>SIMD_EPSILON)
					{
						sepNormalWorldSpace = withoutMargin.m_reportedNormalOnWorld;  //gjkPairDetector.getCachedSeparatingAxis()*(1.f/l2);
						//minDist = -1e30f;//gjkPairDetector.getCachedSeparatingDistance();
						minDist = withoutMargin.m_reportedDistance;  //gjkPairDetector.getCachedSeparatingDistance()+min0->getMargin()+min1->getMargin();

#ifdef ZERO_MARGIN
						foundSepAxis = true;  //gjkPairDetector.getCachedSeparatingDistance()<0.f;
#else
						foundSepAxis = withoutMargin.m_foundResult && minDist < 0;  //-(min0->getMargin()+min1->getMargin());
#endif
					}
				}
				if (foundSepAxis)
				{
					//				printf("sepNormalWorldSpace=%f,%f,%f\n",sepNormalWorldSpace.getX(),sepNormalWorldSpace.getY(),sepNormalWorldSpace.getZ());

					worldVertsB1.resize(0);
					cbtPolyhedralContactClipping::clipHullAgainstHull(sepNormalWorldSpace, *polyhedronA->getConvexPolyhedron(), *polyhedronB->getConvexPolyhedron(),
																	 body0Wrap->getWorldTransform(),
																	 body1Wrap->getWorldTransform(), minDist - threshold, threshold, worldVertsB1, worldVertsB2,
																	 *resultOut);
				}
				if (m_ownManifold)
				{
					resultOut->refreshContactPoints();
				}
				return;
			}
			else
			{
				//we can also deal with convex versus triangle (without connectivity data)
				if (dispatchInfo.m_enableSatConvex && polyhedronA->getConvexPolyhedron() && polyhedronB->getShapeType() == TRIANGLE_SHAPE_PROXYTYPE)
				{
					cbtVertexArray worldSpaceVertices;
					cbtTriangleShape* tri = (cbtTriangleShape*)polyhedronB;
					worldSpaceVertices.push_back(body1Wrap->getWorldTransform() * tri->m_vertices1[0]);
					worldSpaceVertices.push_back(body1Wrap->getWorldTransform() * tri->m_vertices1[1]);
					worldSpaceVertices.push_back(body1Wrap->getWorldTransform() * tri->m_vertices1[2]);

					//tri->initializePolyhedralFeatures();

					cbtScalar threshold = m_manifoldPtr->getContactBreakingThreshold()+ resultOut->m_closestPointDistanceThreshold;

					cbtVector3 sepNormalWorldSpace;
					cbtScalar minDist = -1e30f;
					cbtScalar maxDist = threshold;

					bool foundSepAxis = false;
					bool useSatSepNormal = true;

					if (useSatSepNormal)
					{
#if 0
					if (0)
					{
						//initializePolyhedralFeatures performs a convex hull computation, not needed for a single triangle
						polyhedronB->initializePolyhedralFeatures();
					} else
#endif
						{
							cbtVector3 uniqueEdges[3] = {tri->m_vertices1[1] - tri->m_vertices1[0],
														tri->m_vertices1[2] - tri->m_vertices1[1],
														tri->m_vertices1[0] - tri->m_vertices1[2]};

							uniqueEdges[0].normalize();
							uniqueEdges[1].normalize();
							uniqueEdges[2].normalize();

							cbtConvexPolyhedron polyhedron;
							polyhedron.m_vertices.push_back(tri->m_vertices1[2]);
							polyhedron.m_vertices.push_back(tri->m_vertices1[0]);
							polyhedron.m_vertices.push_back(tri->m_vertices1[1]);

							{
								cbtFace combinedFaceA;
								combinedFaceA.m_indices.push_back(0);
								combinedFaceA.m_indices.push_back(1);
								combinedFaceA.m_indices.push_back(2);
								cbtVector3 faceNormal = uniqueEdges[0].cross(uniqueEdges[1]);
								faceNormal.normalize();
								cbtScalar planeEq = 1e30f;
								for (int v = 0; v < combinedFaceA.m_indices.size(); v++)
								{
									cbtScalar eq = tri->m_vertices1[combinedFaceA.m_indices[v]].dot(faceNormal);
									if (planeEq > eq)
									{
										planeEq = eq;
									}
								}
								combinedFaceA.m_plane[0] = faceNormal[0];
								combinedFaceA.m_plane[1] = faceNormal[1];
								combinedFaceA.m_plane[2] = faceNormal[2];
								combinedFaceA.m_plane[3] = -planeEq;
								polyhedron.m_faces.push_back(combinedFaceA);
							}
							{
								cbtFace combinedFaceB;
								combinedFaceB.m_indices.push_back(0);
								combinedFaceB.m_indices.push_back(2);
								combinedFaceB.m_indices.push_back(1);
								cbtVector3 faceNormal = -uniqueEdges[0].cross(uniqueEdges[1]);
								faceNormal.normalize();
								cbtScalar planeEq = 1e30f;
								for (int v = 0; v < combinedFaceB.m_indices.size(); v++)
								{
									cbtScalar eq = tri->m_vertices1[combinedFaceB.m_indices[v]].dot(faceNormal);
									if (planeEq > eq)
									{
										planeEq = eq;
									}
								}

								combinedFaceB.m_plane[0] = faceNormal[0];
								combinedFaceB.m_plane[1] = faceNormal[1];
								combinedFaceB.m_plane[2] = faceNormal[2];
								combinedFaceB.m_plane[3] = -planeEq;
								polyhedron.m_faces.push_back(combinedFaceB);
							}

							polyhedron.m_uniqueEdges.push_back(uniqueEdges[0]);
							polyhedron.m_uniqueEdges.push_back(uniqueEdges[1]);
							polyhedron.m_uniqueEdges.push_back(uniqueEdges[2]);
							polyhedron.initialize2();

							polyhedronB->setPolyhedralFeatures(polyhedron);
						}

						foundSepAxis = cbtPolyhedralContactClipping::findSeparatingAxis(
							*polyhedronA->getConvexPolyhedron(), *polyhedronB->getConvexPolyhedron(),
							body0Wrap->getWorldTransform(),
							body1Wrap->getWorldTransform(),
							sepNormalWorldSpace, *resultOut);
						//	 printf("sepNormalWorldSpace=%f,%f,%f\n",sepNormalWorldSpace.getX(),sepNormalWorldSpace.getY(),sepNormalWorldSpace.getZ());
					}
					else
					{
#ifdef ZERO_MARGIN
						gjkPairDetector.setIgnoreMargin(true);
						gjkPairDetector.getClosestPoints(input, *resultOut, dispatchInfo.m_debugDraw);
#else
						gjkPairDetector.getClosestPoints(input, dummy, dispatchInfo.m_debugDraw);
#endif  //ZERO_MARGIN

						if (dummy.m_hasContact && dummy.m_depth < 0)
						{
							if (foundSepAxis)
							{
								if (dummy.m_normalOnBInWorld.dot(sepNormalWorldSpace) < 0.99)
								{
									printf("?\n");
								}
							}
							else
							{
								printf("!\n");
							}
							sepNormalWorldSpace.setValue(0, 0, 1);  // = dummy.m_normalOnBInWorld;
							//minDist = dummy.m_depth;
							foundSepAxis = true;
						}
#if 0
					cbtScalar l2 = gjkPairDetector.getCachedSeparatingAxis().length2();
					if (l2>SIMD_EPSILON)
					{
						sepNormalWorldSpace = gjkPairDetector.getCachedSeparatingAxis()*(1.f/l2);
						//minDist = gjkPairDetector.getCachedSeparatingDistance();
						//maxDist = threshold;
						minDist = gjkPairDetector.getCachedSeparatingDistance()-min0->getMargin()-min1->getMargin();
						foundSepAxis = true;
					}
#endif
					}

					if (foundSepAxis)
					{
						worldVertsB2.resize(0);
						cbtPolyhedralContactClipping::clipFaceAgainstHull(sepNormalWorldSpace, *polyhedronA->getConvexPolyhedron(),
																		 body0Wrap->getWorldTransform(), worldSpaceVertices, worldVertsB2, minDist - threshold, maxDist, *resultOut);
					}

					if (m_ownManifold)
					{
						resultOut->refreshContactPoints();
					}

					return;
				}
			}
		}

		/* ***CHRONO*** radial contact normal on the lateral side of cylinder-like shapes */
		const bool lateralPair = cbtLateralNormalResult::applies(min0, min1);
		cbtLateralNormalResult lateralResult(resultOut, min0, min1, input.m_transformA, input.m_transformB);
		if (lateralPair)
			gjkPairDetector.getClosestPoints(input, lateralResult, dispatchInfo.m_debugDraw);
		else
			gjkPairDetector.getClosestPoints(input, *resultOut, dispatchInfo.m_debugDraw);

		//now perform 'm_numPerturbationIterations' collision queries with the perturbated collision objects

		//perform perturbation when more then 'm_minimumPointsPerturbationThreshold' points
		/* ***CHRONO*** not for a strictly convex shape (sphere, ellipsoid), which touches in a single point, nor when the
		shape to perturb (the one with the smaller angular motion disc, see below) is an element of a mesh, beam or node
		cloud */
		if (m_numPerturbationIterations && !strictlyConvexPair && resultOut->getPersistentManifold()->getNumContacts() < m_minimumPointsPerturbationThreshold &&
			!isMeshElement(min0->getAngularMotionDisc() < min1->getAngularMotionDisc() ? min0 : min1))
		{
			int i;
			cbtVector3 v0, v1;
			cbtVector3 sepNormalWorldSpace;
			cbtScalar l2 = gjkPairDetector.getCachedSeparatingAxis().length2();

			if (l2 > SIMD_EPSILON)
			{
				sepNormalWorldSpace = gjkPairDetector.getCachedSeparatingAxis() * (1.f / l2);
				/* ***CHRONO*** perturbed points are compared with the corrected normal */
				cbtVector3 filterNormal = lateralResult.m_corrected ? lateralResult.m_normal : sepNormalWorldSpace.normalized();

				cbtPlaneSpace1(sepNormalWorldSpace, v0, v1);

				bool perturbeA = true;
				const cbtScalar angleLimit = 0.125f * SIMD_PI;
				cbtScalar perturbeAngle;
				cbtScalar radiusA = min0->getAngularMotionDisc();
				cbtScalar radiusB = min1->getAngularMotionDisc();
				if (radiusA < radiusB)
				{
					perturbeAngle = gContactBreakingThreshold / radiusA;
					perturbeA = true;
				}
				else
				{
					perturbeAngle = gContactBreakingThreshold / radiusB;
					perturbeA = false;
				}
				if (perturbeAngle > angleLimit)
					perturbeAngle = angleLimit;

				cbtTransform unPerturbedTransform;
				if (perturbeA)
				{
					unPerturbedTransform = input.m_transformA;
				}
				else
				{
					unPerturbedTransform = input.m_transformB;
				}

				for (i = 0; i < m_numPerturbationIterations; i++)
				{
					if (v0.length2() > SIMD_EPSILON)
					{
						cbtQuaternion perturbeRot(v0, perturbeAngle);
						cbtScalar iterationAngle = i * (SIMD_2_PI / cbtScalar(m_numPerturbationIterations));
						cbtQuaternion rotq(sepNormalWorldSpace, iterationAngle);

						if (perturbeA)
						{
							input.m_transformA.setBasis(cbtMatrix3x3(rotq.inverse() * perturbeRot * rotq) * body0Wrap->getWorldTransform().getBasis());
							input.m_transformB = body1Wrap->getWorldTransform();
#ifdef DEBUG_CONTACTS
							dispatchInfo.m_debugDraw->drawTransform(input.m_transformA, 10.0);
#endif  //DEBUG_CONTACTS
						}
						else
						{
							input.m_transformA = body0Wrap->getWorldTransform();
							input.m_transformB.setBasis(cbtMatrix3x3(rotq.inverse() * perturbeRot * rotq) * body1Wrap->getWorldTransform().getBasis());
#ifdef DEBUG_CONTACTS
							dispatchInfo.m_debugDraw->drawTransform(input.m_transformB, 10.0);
#endif
						}

						/* ***CHRONO*** normals within a tenth of the perturbation angle of the unperturbed one (but not below the
						rounding errors of the normals) */
						cbtScalar normalTolerance = cbtMax(cbtScalar(1) - cbtCos(cbtScalar(0.1) * perturbeAngle), cbtScalar(8) * SIMD_EPSILON);
						cbtPerturbedContactResult perturbedResultOut(resultOut, input.m_transformA, input.m_transformB, unPerturbedTransform, perturbeA, dispatchInfo.m_debugDraw,
																	perturbeA ? min0->getMargin() : min1->getMargin(), filterNormal, normalTolerance);
						if (lateralPair)
						{
							cbtLateralNormalResult lateralPerturbed(&perturbedResultOut, min0, min1, input.m_transformA, input.m_transformB);
							gjkPairDetector.getClosestPoints(input, lateralPerturbed, dispatchInfo.m_debugDraw);
						}
						else
							gjkPairDetector.getClosestPoints(input, perturbedResultOut, dispatchInfo.m_debugDraw);
					}
				}
			}
		}

#ifdef USE_SEPDISTANCE_UTIL2
		if (dispatchInfo.m_useConvexConservativeDistanceUtil && (sepDist > SIMD_EPSILON))
		{
			m_sepDistance.initSeparatingDistance(gjkPairDetector.getCachedSeparatingAxis(), sepDist, body0->getWorldTransform(), body1->getWorldTransform());
		}
#endif  //USE_SEPDISTANCE_UTIL2
	}

	if (m_ownManifold)
	{
		resultOut->refreshContactPoints();
	}
}

bool disableCcd = false;
cbtScalar cbtConvexConvexAlgorithm::calculateTimeOfImpact(cbtCollisionObject* col0, cbtCollisionObject* col1, const cbtDispatcherInfo& dispatchInfo, cbtManifoldResult* resultOut)
{
	(void)resultOut;
	(void)dispatchInfo;
	///Rather then checking ALL pairs, only calculate TOI when motion exceeds threshold

	///Linear motion for one of objects needs to exceed m_ccdSquareMotionThreshold
	///col0->m_worldTransform,
	cbtScalar resultFraction = cbtScalar(1.);

	cbtScalar squareMot0 = (col0->getInterpolationWorldTransform().getOrigin() - col0->getWorldTransform().getOrigin()).length2();
	cbtScalar squareMot1 = (col1->getInterpolationWorldTransform().getOrigin() - col1->getWorldTransform().getOrigin()).length2();

	if (squareMot0 < col0->getCcdSquareMotionThreshold() &&
		squareMot1 < col1->getCcdSquareMotionThreshold())
		return resultFraction;

	if (disableCcd)
		return cbtScalar(1.);

	//An adhoc way of testing the Continuous Collision Detection algorithms
	//One object is approximated as a sphere, to simplify things
	//Starting in penetration should report no time of impact
	//For proper CCD, better accuracy and handling of 'allowed' penetration should be added
	//also the mainloop of the physics should have a kind of toi queue (something like Brian Mirtich's application of Timewarp for Rigidbodies)

	/// Convex0 against sphere for Convex1
	{
		cbtConvexShape* convex0 = static_cast<cbtConvexShape*>(col0->getCollisionShape());

		cbtSphereShape sphere1(col1->getCcdSweptSphereRadius());  //todo: allow non-zero sphere sizes, for better approximation
		cbtConvexCast::CastResult result;
		cbtVoronoiSimplexSolver voronoiSimplex;
		//SubsimplexConvexCast ccd0(&sphere,min0,&voronoiSimplex);
		///Simplification, one object is simplified as a sphere
		cbtGjkConvexCast ccd1(convex0, &sphere1, &voronoiSimplex);
		//ContinuousConvexCollision ccd(min0,min1,&voronoiSimplex,0);
		if (ccd1.calcTimeOfImpact(col0->getWorldTransform(), col0->getInterpolationWorldTransform(),
								  col1->getWorldTransform(), col1->getInterpolationWorldTransform(), result))
		{
			//store result.m_fraction in both bodies

			if (col0->getHitFraction() > result.m_fraction)
				col0->setHitFraction(result.m_fraction);

			if (col1->getHitFraction() > result.m_fraction)
				col1->setHitFraction(result.m_fraction);

			if (resultFraction > result.m_fraction)
				resultFraction = result.m_fraction;
		}
	}

	/// Sphere (for convex0) against Convex1
	{
		cbtConvexShape* convex1 = static_cast<cbtConvexShape*>(col1->getCollisionShape());

		cbtSphereShape sphere0(col0->getCcdSweptSphereRadius());  //todo: allow non-zero sphere sizes, for better approximation
		cbtConvexCast::CastResult result;
		cbtVoronoiSimplexSolver voronoiSimplex;
		//SubsimplexConvexCast ccd0(&sphere,min0,&voronoiSimplex);
		///Simplification, one object is simplified as a sphere
		cbtGjkConvexCast ccd1(&sphere0, convex1, &voronoiSimplex);
		//ContinuousConvexCollision ccd(min0,min1,&voronoiSimplex,0);
		if (ccd1.calcTimeOfImpact(col0->getWorldTransform(), col0->getInterpolationWorldTransform(),
								  col1->getWorldTransform(), col1->getInterpolationWorldTransform(), result))
		{
			//store result.m_fraction in both bodies

			if (col0->getHitFraction() > result.m_fraction)
				col0->setHitFraction(result.m_fraction);

			if (col1->getHitFraction() > result.m_fraction)
				col1->setHitFraction(result.m_fraction);

			if (resultFraction > result.m_fraction)
				resultFraction = result.m_fraction;
		}
	}

	return resultFraction;
}
