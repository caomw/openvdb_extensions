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


#include <NearestNeighborOperation.h>
using namespace openvdb;
void NearestNeighborOperation::init(SpringLevelSet& mGrid) {
	NearestNeighborMap& map = mGrid.mNearestNeighbors;
	map.clear();
	map.resize(mGrid.mConstellation.getNumVertexes(),
			std::list<SpringlNeighbor>());
}
void NearestNeighborOperation::compute(Springl& springl, SpringLevelSet& mGrid,
		double t) {
	const float D2 = SpringLevelSet::NEAREST_NEIGHBOR_RANGE
			* SpringLevelSet::NEAREST_NEIGHBOR_RANGE;

	openvdb::math::DenseStencil<openvdb::Int32Grid> stencil =
			openvdb::math::DenseStencil<openvdb::Int32Grid>(
					*mGrid.mSpringlIndexGrid,
					ceil(SpringLevelSet::NEAREST_NEIGHBOR_RANGE));
	openvdb::Vec3s refPoint = springl.particle();

	stencil.moveTo(
			Coord(std::floor(refPoint[0] + 0.5f),
					std::floor(refPoint[1] + 0.5f),
					std::floor(refPoint[2] + 0.5f)));
	int sz = stencil.size();
	if (sz == 0)
		return;
	Index32 N = mGrid.mConstellation.getNumSpringls();
	std::vector<openvdb::Index32> stencilCopy;
	for (int i = 0; i < sz; i++) {
		openvdb::Index32 id = stencil.getValue(i);
		if (id >= N)
			continue;
		if (id != springl.id) {
			stencilCopy.push_back(id);
		}
	}
	if (stencilCopy.size() == 0)
		return;
	std::sort(stencilCopy.begin(), stencilCopy.end());
	sz = stencilCopy.size();

	openvdb::Index32 last = -1;
	SpringlNeighbor bestNbr;
	std::vector<SpringlNeighbor> tmpRange;
	for (int k = 0; k < springl.size(); k++) {
		std::list<SpringlNeighbor>& mapList = mGrid.getNearestNeighbors(
				springl.id, k);
		refPoint = springl[k];
		//
		last = -1;

		tmpRange.clear();
		for (int i = 0; i < sz; i++) {
			openvdb::Index32 nbrId = stencilCopy[i];
			if (nbrId == last)
				continue;
			Springl& snbr = mGrid.getSpringl(nbrId);
			bestNbr = SpringlNeighbor(nbrId, -1, D2);
			for (int8_t n = 0; n < snbr.size(); n++) {
				float d = snbr.distanceToEdgeSqr(refPoint, n);
				if (d <= bestNbr.distance) {
					bestNbr.edgeId = n;
					bestNbr.distance = d;
				}
			}
			if (bestNbr.edgeId >= 0)
				tmpRange.push_back(bestNbr);
			last = nbrId;
		}
		sort(tmpRange.begin(), tmpRange.end());
		for (int nn = 0, nmax = std::min(SpringLevelSet::MAX_NEAREST_NEIGHBORS,
				(int) tmpRange.size()); nn < nmax; nn++) {
			mapList.push_back(tmpRange[nn]);
		}
		//if (bestNbr.springlId >= 0 && bestNbr.springlId < N)mapList.push_back(bestNbr);
	}
}
