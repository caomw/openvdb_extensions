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

#include "Springl.h"
#include "Constellation.h"
#include "GeometryUtil.h"
#include <AlloyMath.h>
#include <openvdb/util/Util.h>
using namespace openvdb;
openvdb::Vec3s& Springl::normal() const {
	return mesh->mParticleNormals[id];
}
openvdb::Vec3s& Springl::particle() const {
	return mesh->mParticles[id];
}
openvdb::Vec3s& Springl::particleVelocity() const {
	return mesh->mParticleVelocity[id];
}
openvdb::Vec3s& Springl::vertexVelocity(size_t idx) const {
	return mesh->mVertexVelocity[offset + idx];
}
uint8_t& Springl::label() const {
	return mesh->mParticleLabel[id];
}
openvdb::Vec3s& Springl::operator[](size_t idx) {
	return mesh->mVertexes[offset + idx];
}
const openvdb::Vec3s& Springl::operator[](size_t idx) const {
	return mesh->mVertexes[offset + idx];
}
Springl::Springl(Constellation* _mesh) :
		 mesh(_mesh) ,id(0), offset(0){
}
std::ostream& operator<<(std::ostream& ostr, const SpringlNeighbor& classname) {
	ostr << "{" << classname.springlId << "|" << static_cast<int>(classname.edgeId) << ":" << std::setprecision(4) << classname.distance << "}";
	return ostr;
}
int Springl::size() const {
	return (id<mesh->mQuads.size())?4:3;//This works because quad indexes are assigned first in springl list.
}

float Springl::distanceToParticle(const openvdb::Vec3s& pt) {
	return ((particle()) - pt).length();
}

float Springl::distanceToParticleSqr(const openvdb::Vec3s& pt) {
	return ((particle()) - pt).lengthSqr();
}
float Springl::distanceToFace(const openvdb::Vec3s& pt) {
	return std::sqrt(distanceToFaceSqr(pt));
}

float Springl::distanceToFaceSqr(const openvdb::Vec3s& pt) {
	Vec3s closest;
	if (size() == 3) {
		return DistanceToTriangleSqr(pt, (*this)[0], (*this)[1], (*this)[2], &closest);
	} else {
		return DistanceToQuadSqr(pt, (*this)[0], (*this)[1], (*this)[2], (*this)[3], normal(), &closest);
	}
}
float Springl::signedDistanceToFaceSqr(const openvdb::Vec3s& pt) {
	Vec3s closest;
	if (size() == 3) {
		float d = DistanceToTriangleSqr(pt, (*this)[0], (*this)[1], (*this)[2], &closest);
		return d * aly::sign((pt - closest).dot(normal()));

	} else {
		float d = DistanceToQuadSqr(pt, (*this)[0], (*this)[1], (*this)[2], (*this)[3], normal(), &closest);
		return d * aly::sign((pt - closest).dot(normal()));
	}
}
openvdb::math::BBox<Vec3s> Springl::getBoundingBox() {
	openvdb::math::BBox<Vec3s> bbox;

	bbox.min() = Vec3s(std::numeric_limits<float>::max());
	bbox.max() = Vec3s(std::numeric_limits<float>::min());

	for (int n = 0; n < size(); n++) {
		Vec3s pt = (*this)[n];
		bbox.min() = openvdb::math::Min(bbox.min(), pt);
		bbox.max() = openvdb::math::Max(bbox.max(), pt);
	}
	return bbox;
}
float Springl::signedDistanceToFace(const openvdb::Vec3s& pt) {
	float d = signedDistanceToFaceSqr(pt);
	return aly::sign(d) * std::sqrt(abs(d));
}
float Springl::distanceToEdgeSqr(const openvdb::Vec3s& pt, int e) {
	return DistanceToEdgeSqr(pt, (*this)[e], (*this)[(e + 1) % size()]);
}
float Springl::distanceToEdge(const openvdb::Vec3s& pt, int e) {
	return std::sqrt(DistanceToEdgeSqr(pt, (*this)[e], (*this)[(e + 1) % size()]));
}

openvdb::Vec3s Springl::computeCentroid() const {
	openvdb::Vec3s centroid = openvdb::Vec3s(0.0f, 0.0f, 0.0f);
	int K = size();
	for (int k = 0; k < K; k++) {
		centroid += (*this)[k];
	}
	centroid = (1.0 / K) * centroid;
	return centroid;
}
openvdb::Vec3s Springl::computeNormal(const float eps) const {
	openvdb::Vec3s norm(0.0f);
	int K = size();
	openvdb::Vec3s pt = particle();
	for (int k = 0; k < K; k++) {
		norm += ((*this)[k] - pt).cross((*this)[(k + 1) % K] - pt);
	}

	norm.normalize(eps);
	return norm;
}
float Springl::area() const {
	openvdb::Vec3s norm(0.0f);
	int K = size();
	for (int k = 0; k < K; k++) {
		norm += ((*this)[k] - (particle())).cross((*this)[(k + 1) % K] - (particle()));
	}
	return 0.5f * norm.length();
}
