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

#include "RelaxOperation.h"
#include "GeometryUtil.h"
using namespace openvdb;
void RelaxOperation::init(SpringLevelSet& mGrid) {
	mGrid.mConstellation.mVertexAuxBuffer.resize(mGrid.mConstellation.getNumVertexes());
}
void RelaxOperation::apply(Springl& springl, SpringLevelSet& mGrid, double dt) {
	int K = springl.size();
	for (int k = 0; k < K; k++) {
		springl[k] = mGrid.mConstellation.mVertexAuxBuffer[springl.offset + k];
	}
}
void RelaxOperation::compute(Springl& springl, SpringLevelSet& mGrid, double t) {
	float w, len;
	Vec3s tanget;
	Vec3s dir;
	int K = springl.size();
	std::vector<Vec3s> vertexVelocity(K);
	std::vector<Vec3s> tangets(K);
	std::vector<float> springForce(K);
	std::vector<float> tangetLengths(K);
	Vec3s particlePt = springl.particle();
	Vec3s startVelocity = Vec3s(0);
	Vec3s resultantMoment = Vec3s(0);
	const float MAX_FORCE = 0.999f;
	Vec3s start;
	float dotProd;
	Vec3s pt2;
	for (int k = 0; k < K; k++) {
		std::list<SpringlNeighbor>& map = mGrid.getNearestNeighbors(springl.id, k);
		start = springl[k];
		// edge from pivot to magnet
		tanget = (start - particlePt);
		tangetLengths[k] = tanget.length();
		if (tangetLengths[k] > 1E-6f) {
			tanget *= (1.0f / tangetLengths[k]);
		}
		tangets[k] = tanget;
		startVelocity = Vec3s(0);
		// Sum forces
		//unroll loop
		for (SpringlNeighbor ci : map) {
			//Closest point should be recomputed each time and does not need to be stored

			Springl& nbr = mGrid.getSpringl(ci.springlId);
			DistanceToEdgeSqr(start, nbr[ci.edgeId], nbr[(ci.edgeId + 1) % nbr.size()], &pt2);
			dir = (pt2 - start);
			len = dir.length();
			w = ((len - 2 * SpringLevelSet::PARTICLE_RADIUS) / (SpringLevelSet::MAX_VEXT + 2 * SpringLevelSet::PARTICLE_RADIUS));
			w = atanh(MAX_FORCE * aly::clamp(w, -1.0f, 1.0f));
			startVelocity += (w * dir);
		}
		if (map.size() > 0)
			startVelocity /= map.size();
		vertexVelocity[k] = SpringLevelSet::RELAX_TIMESTEP * startVelocity * SpringLevelSet::SHARPNESS;
		springForce[k] = SpringLevelSet::RELAX_TIMESTEP * SpringLevelSet::SPRING_CONSTANT * (2 * SpringLevelSet::PARTICLE_RADIUS - tangetLengths[k]);
		resultantMoment += vertexVelocity[k].cross(tangets[k]);
	}
	//std::cout<<"moment "<<resultantMoment<<" Normal "<<*springl.normal<<std::endl;

	openvdb::math::Mat3<float> rot = CreateAxisAngle(resultantMoment, -resultantMoment.length());

	//std::cout<<"Rotation\n"<<rot<<std::endl;

	//std::cout<<springl.id<<springl.offset<<" ROTATION "<<resultantMoment<<" "<<springForce[0]<<" "<<vertexVelocity[0]<<std::endl;
	for (int k = 0; k < K; k++) {
		start = springl[k] - particlePt;
		dotProd = std::max(start.length() + vertexVelocity[k].dot(tangets[k]) + springForce[k], 0.001f);
		start = dotProd * tangets[k];

		//disable rotation
		start = rot * start;
		mGrid.mConstellation.mVertexAuxBuffer[springl.offset + k] = start + particlePt;
	}
}

