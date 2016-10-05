/*
 * Copyright(C) 2016, Blake C. Lucas, Ph.D. (img.science@gmail.com)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
#include "GeometryUtil.h"
using namespace openvdb;
openvdb::math::Mat3<float> CreateAxisAngle(Vec3s a1, float angle) {
	openvdb::math::Mat3<float> M;
	float mag = a1.length();
	if (mag < 1E-6f) {
		M[0][0] = 1.0f;
		M[0][1] = 0.0f;
		M[0][2] = 0.0f;

		M[1][0] = 0.0f;
		M[1][1] = 1.0f;
		M[1][2] = 0.0f;

		M[2][0] = 0.0f;
		M[2][1] = 0.0f;
		M[2][2] = 1.0f;
	} else {
		mag = 1.0f / mag;
		float ax = a1[0] * mag;
		float ay = a1[1] * mag;
		float az = a1[2] * mag;
		float sinTheta = (float) sin(angle);
		float cosTheta = (float) cos(angle);
		float t = 1.0f - cosTheta;

		float xz = ax * az;
		float xy = ax * ay;
		float yz = ay * az;

		M[0][0] = t * ax * ax + cosTheta;
		M[0][1] = t * xy - sinTheta * az;
		M[0][2] = t * xz + sinTheta * ay;

		M[1][0] = t * xy + sinTheta * az;
		M[1][1] = t * ay * ay + cosTheta;
		M[1][2] = t * yz - sinTheta * ax;

		M[2][0] = t * xz - sinTheta * ay;
		M[2][1] = t * yz + sinTheta * ax;
		M[2][2] = t * az * az + cosTheta;
	}
	return M;
}
float DistanceToEdgeSqr(const openvdb::Vec3s& pt, const openvdb::Vec3s& pt1, const openvdb::Vec3s& pt2, openvdb::Vec3s* lastClosestSegmentPoint) {
	using namespace openvdb;
	using namespace openvdb::math;
	Vec3s dir = pt2 - pt1;
	float len = dir.length();
	dir.normalize(1E-6f);
	Vec3s diff = pt - pt1;
	float mSegmentParameter = dir.dot(diff);
	if (0 < mSegmentParameter) {
		if (mSegmentParameter < len) {
			*lastClosestSegmentPoint = dir * mSegmentParameter + pt1;
		} else {
			*lastClosestSegmentPoint = pt2;
		}
	} else {
		*lastClosestSegmentPoint = pt1;
	}
	return (pt - (*lastClosestSegmentPoint)).lengthSqr();
}
float DistanceToEdgeSqr(const openvdb::Vec3s& pt, const openvdb::Vec3s& pt1, const openvdb::Vec3s& pt2) {
	openvdb::Vec3s tmp;
	return DistanceToEdgeSqr(pt, pt1, pt2, &tmp);
}
//Implementation from geometric tools (http://www.geometrictools.com)
inline Vec3s parametricTriangle(Vec3s e0, Vec3s e1, float s, float t, Vec3s B) {
	Vec3s Bsum = B + s * e0 + t * e1;
	return Bsum;
}
float Angle(openvdb::Vec3s& v0, openvdb::Vec3s& v1, openvdb::Vec3s& v2) {
	Vec3s v = v0 - v1;
	Vec3s w = v2 - v1;
	float len1 = v.length();
	float len2 = w.length();
	return std::acos(v.dot(w) / std::max(1E-8f, len1 * len2));
}
float DistanceToTriangleSqr(const openvdb::Vec3s& p, const openvdb::Vec3s& v0, const openvdb::Vec3s& v1, const openvdb::Vec3s& v2,
		openvdb::Vec3s* closestPoint) {
	float distanceSquared = 0;
	int region_id = 0;

	Vec3s P = p;
	Vec3s B = v0;
	Vec3s e0 = v1 - v0;
	Vec3s e1 = v2 - v0;
	float a = e0.dot(e0);
	float b = e0.dot(e1);
	float c = e1.dot(e1);
	Vec3s dv = B - P;
	float d = e0.dot(dv);
	float e = e1.dot(dv);
	// Determine which region_id contains s, t

	float det = a * c - b * b;
	float s = b * e - c * d;
	float t = b * d - a * e;

	if (s + t <= det) {
		if (s < 0) {
			if (t < 0) {
				region_id = 4;
			} else {
				region_id = 3;
			}
		} else if (t < 0) {
			region_id = 5;
		} else {
			region_id = 0;
		}
	} else {
		if (s < 0) {
			region_id = 2;
		} else if (t < 0) {
			region_id = 6;
		} else {
			region_id = 1;
		}
	}

	// Parametric Triangle Point
	Vec3s T(0.0f);

	if (region_id == 0) { // Region 0
		float invDet = (float) 1 / (float) det;
		s *= invDet;
		t *= invDet;

		// Find point on parametric triangle based on s and t
		T = parametricTriangle(e0, e1, s, t, B);
		// Find distance from P to T
		Vec3s tmp = P - T;
		distanceSquared = tmp.lengthSqr();
	} else if (region_id == 1) { // Region 1
		float numer = c + e - b - d;

		if (numer < +0) {
			s = 0;
		} else {
			float denom = a - 2 * b + c;
			s = (numer >= denom ? 1 : numer / denom);
		}
		t = 1 - s;

		// Find point on parametric triangle based on s and t
		T = parametricTriangle(e0, e1, s, t, B);
		// Find distance from P to T
		Vec3s tmp = P - T;
		distanceSquared = tmp.lengthSqr();
	} else if (region_id == 2) { // Region 2
		float tmp0 = b + d;
		float tmp1 = c + e;

		if (tmp1 > tmp0) {
			float numer = tmp1 - tmp0;
			float denom = a - 2 * b + c;
			s = (numer >= denom ? 1 : numer / denom);
			t = 1 - s;
		} else {
			s = 0;
			t = (tmp1 <= 0 ? 1 : (e >= 0 ? 0 : -e / c));
		}

		// Find point on parametric triangle based on s and t
		T = parametricTriangle(e0, e1, s, t, B);
		// Find distance from P to T
		Vec3s tmp = P - T;
		distanceSquared = tmp.lengthSqr();
	} else if (region_id == 3) { // Region 3
		s = 0;
		t = (e >= 0 ? 0 : (-e >= c ? 1 : -e / c));

		// Find point on parametric triangle based on s and t
		T = parametricTriangle(e0, e1, s, t, B);
		// Find distance from P to T
		Vec3s tmp = P - T;
		distanceSquared = tmp.lengthSqr();

	} else if (region_id == 4) { // Region 4
		float tmp0 = c + e;
		float tmp1 = a + d;

		if (tmp0 > tmp1) {
			s = 0;
			t = (tmp1 <= 0 ? 1 : (e >= 0 ? 0 : -e / c));
		} else {
			t = 0;
			s = (tmp1 <= 0 ? 1 : (d >= 0 ? 0 : -d / a));
		}

		// Find point on parametric triangle based on s and t
		T = parametricTriangle(e0, e1, s, t, B);
		// Find distance from P to T
		Vec3s tmp = P - T;
		distanceSquared = tmp.lengthSqr();
	} else if (region_id == 5) { // Region 5
		t = 0;
		s = (d >= 0 ? 0 : (-d >= a ? 1 : -d / a));

		// Find point on parametric triangle based on s and t
		T = parametricTriangle(e0, e1, s, t, B);
		// Find distance from P to T
		Vec3s tmp = P - T;
		distanceSquared = tmp.lengthSqr();
	} else { // Region 6
		float tmp0 = b + e;
		float tmp1 = a + d;

		if (tmp1 > tmp0) {
			float numer = tmp1 - tmp0;
			float denom = c - 2 * b + a;
			t = (numer >= denom ? 1 : numer / denom);
			s = 1 - t;
		} else {
			t = 0;
			s = (tmp1 <= 0 ? 1 : (d >= 0 ? 0 : -d / a));
		}

		// Find point on parametric triangle based on s and t
		T = parametricTriangle(e0, e1, s, t, B);
		// Find distance from P to T
		Vec3s tmp = P - T;
		distanceSquared = tmp.lengthSqr();
	}
	(*closestPoint) = T;
	return distanceSquared;
}

//What if quad is non-convex? Does this hold?
float DistanceToQuadSqr(const openvdb::Vec3s& p, const openvdb::Vec3s& v0, const openvdb::Vec3s& v1, const openvdb::Vec3s& v2, const openvdb::Vec3s& v3,
		const openvdb::Vec3s& norm, openvdb::Vec3s* closestPoint) {
	Vec3s cp1;
	Vec3s cp2;
	float d1, d2;
	if ((v2 - v0).cross(v1 - v0).dot(norm) > 0) {
		d1 = DistanceToTriangleSqr(p, v0, v1, v2, &cp1);
		d2 = DistanceToTriangleSqr(p, v2, v3, v0, &cp2);
	} else {
		d1 = DistanceToTriangleSqr(p, v1, v2, v3, &cp1);
		d2 = DistanceToTriangleSqr(p, v3, v0, v1, &cp2);
	}

	if (d1 < d2) {
		*closestPoint = cp1;
		return d1;
	} else {
		*closestPoint = cp2;
		return d2;
	}
}
