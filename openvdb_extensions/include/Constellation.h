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
#ifndef INCLUDE_CONSTELLATION_H_
#define INCLUDE_CONSTELLATION_H_
#include "Springl.h"
#include <openvdb/openvdb.h>
namespace aly{
	class Mesh;
}
class Constellation {
protected:
	openvdb::math::Mat4f mPose;
	openvdb::math::BBox<openvdb::Vec3d> mBoundingBox;
public:
	std::vector<Springl> springls;
	float mMaxVelocityMagnitude = 0;
	float mMinVelocityMagnitude = 0;
	std::vector<openvdb::Vec3s> mLines;
	std::vector<openvdb::Vec3s> mParticles;
	std::vector<openvdb::Vec3s> mParticleNormals;
	std::vector<openvdb::Vec3s> mVertexes;
	std::vector<openvdb::Vec3s> mColors;
	std::vector<openvdb::Vec3s> mVertexNormals;
	std::vector<openvdb::Vec2s> textureMap;
	std::vector<openvdb::Vec3s> mVertexAuxBuffer;
	std::vector<openvdb::Vec3s> mParticleVelocity;
	std::vector<openvdb::Vec3s> mVertexVelocity;
	std::vector<uint8_t> mParticleLabel;
	std::vector<openvdb::Vec4I> mQuads;
	std::vector<openvdb::Vec3I> mTriangles;
	Constellation();
	void initialize();
	inline openvdb::BBoxd getBoundingBox() {
		return mBoundingBox;
	}
	openvdb::math::BBox<openvdb::Vec3d>& updateBoundingBox();
	void scale(float sc);
	inline void setPose(openvdb::Mat4s& pose) {
		mPose = pose;
	}
	inline openvdb::Mat4s& getPose() {
		return mPose;
	}
	inline size_t getNumSpringls() const {
		return springls.size();
	}
	inline size_t getNumVertexes() const {
		return mVertexes.size();
	}

	Springl& operator[](size_t idx) {
		return springls[idx];
	}
	const Springl& operator[](size_t idx) const {
		return springls[idx];
	}
	openvdb::Vec3s closestPointOnEdge(const openvdb::Vec3s& start,
			const SpringlNeighbor& ci);
	void dilate(float distance);
	void reset();
	void copyFrom(const Constellation& mesh);
	void copyTo(aly::Mesh& mesh);
	void updateVertexNormals(int SMOOTH_ITERATIONS = 0, float DOT_TOLERANCE = 0.75f);
	void mapIntoBoundingBox(float voxelSize);
	void mapOutOfBoundingBox(float voxelSize);
	float estimateVoxelSize(int stride = 4);
	~Constellation();
};
bool WriteConstellationToFile(const std::string& file,const Constellation& mesh);
bool ReadConstellationFromFile(const std::string& file,Constellation& mesh);
#endif /* INCLUDE_CONSTELLATION_H_ */
