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

#ifndef INCLUDE_SPRINGL_H_
#define INCLUDE_SPRINGL_H_
#include <openvdb/openvdb.h>
class Constellation;
class Springl {
private:
	Constellation* mesh;
public:
	openvdb::Index32 id;
	openvdb::Index32 offset;
	openvdb::Vec3s& normal() const;
	openvdb::Vec3s& particle() const;
	openvdb::Vec3s& particleVelocity() const;
	openvdb::Vec3s& vertexVelocity(size_t idx) const;
	uint8_t& label() const;
	openvdb::Vec3s& operator[](size_t idx);
	const openvdb::Vec3s& operator[](size_t idx) const ;
	Springl(Constellation* _mesh = nullptr) ;
	int size() const;
	float area() const;
	openvdb::math::BBox<openvdb::Vec3s> getBoundingBox();
	float distanceToFace(const openvdb::Vec3s& pt);
	float signedDistanceToFace(const openvdb::Vec3s& pt);
	float distanceToFaceSqr(const openvdb::Vec3s& pt);
	float signedDistanceToFaceSqr(const openvdb::Vec3s& pt);
	float distanceToParticle(const openvdb::Vec3s& pt);
	float distanceToParticleSqr(const openvdb::Vec3s& pt);
	float distanceToEdgeSqr(const openvdb::Vec3s& pt, int e);
	float distanceToEdge(const openvdb::Vec3s& pt, int e);
	openvdb::Vec3s computeCentroid() const;
	openvdb::Vec3s computeNormal(const float eps = 1E-6f) const;
	~Springl() {
	}
};
struct SpringlNeighbor {
public:
	openvdb::Index32 springlId;
	int edgeId;
	float distance;
	SpringlNeighbor(openvdb::Index32 id = 0, int nbr = -1, float _distance = 0) :
			springlId(id), edgeId(nbr), distance(_distance) {
	}

	friend bool operator<(const SpringlNeighbor& first,
			const SpringlNeighbor& second) {
		return (first.distance < second.distance);
	}
};
std::ostream& operator<<(std::ostream& ostr, const SpringlNeighbor& classname);
#endif /* INCLUDE_SPRINGL_H_ */
