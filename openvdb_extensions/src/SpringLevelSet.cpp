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
#include "SpringLevelSet.h"
#include "NearestNeighborOperation.h"
#include "RelaxOperation.h"
#include "AdvectionForce.h"
#include "UtilitiesVDB.h"
#include <openvdb/Grid.h>
#include <openvdb/tree/LeafManager.h>
#include <openvdb/tools/ChangeBackground.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/LevelSetAdvect.h>
#include <fstream>
#include <cereal/archives/json.hpp>
using namespace openvdb;
typedef openvdb::tools::DiscreteField<openvdb::VectorGrid> VelocityField;
typedef openvdb::tools::LevelSetAdvection<openvdb::FloatGrid, VelocityField> AdvectionTool;
const float SpringLevelSet::NEAREST_NEIGHBOR_RANGE = 1.5f;
const float SpringLevelSet::PARTICLE_RADIUS = 0.05f;
const float SpringLevelSet::MAX_VEXT = 0.5f;
const int SpringLevelSet::MAX_NEAREST_NEIGHBORS = 2;
const float SpringLevelSet::FILL_DISTANCE = 0.3f;
const float SpringLevelSet::CLEAN_DISTANCE = 0.625f;
const float SpringLevelSet::SHARPNESS = 5.0f;
const float SpringLevelSet::SPRING_CONSTANT = 0.3f;
const float SpringLevelSet::RELAX_TIMESTEP = 0.1f;
const float SpringLevelSet::MIN_AREA = 0.05f;
const float SpringLevelSet::MAX_AREA = 2.0;
const float SpringLevelSet::MIN_ASPECT_RATIO = 0.1f;
void SpringLevelSet::updateNearestNeighbors(bool threaded) {
	using namespace openvdb;
	::NearestNeighbors<openvdb::util::NullInterrupter> nn(*this);
	nn.process();
}

openvdb::Vec3s& SpringLevelSet::getParticle(const openvdb::Index32 id) {
	return (mConstellation.springls[id].particle());
}
openvdb::Vec3s& SpringLevelSet::getParticleNormal(const openvdb::Index32 id) {
	return (mConstellation.springls[id].normal());
}
openvdb::Vec3s& SpringLevelSet::getSpringlVertex(const openvdb::Index32 id, const int i) {
	return mConstellation.springls[id][i];
}
openvdb::Vec3s& SpringLevelSet::getSpringlVertex(const openvdb::Index32 id) {
	return mConstellation.mVertexes[id];
}
Springl& SpringLevelSet::getSpringl(const openvdb::Index32 id) {
	return mConstellation.springls[id];
}
void SpringLevelSet::updateLines() {
	std::vector<Vec3s>& lines = mConstellation.mLines;
	lines.clear();
	Vec3s pt, qt;
	for (Index32 i = 0; i < mConstellation.getNumSpringls(); i++) {

		Springl& springl = mConstellation.springls[i];

		for (int k = 0; k < springl.size(); k++) {
			for (SpringlNeighbor nbr : getNearestNeighbors(i, k)) {
				pt = springl[k];
				lines.push_back(pt);
				qt = mConstellation.closestPointOnEdge(pt, nbr);
				lines.push_back(qt);
			}
		}
	}
}
void SpringLevelSet::write(const std::string& file) {
	saveResources(aly::RemoveTrailingSlash(aly::GetParentDirectory(file)));
	std::ofstream os(file);
	cereal::JSONOutputArchive archive(os);
	archive(cereal::make_nvp("springls", *this));
}
void SpringLevelSet::read(const std::string& file) {
	std::ifstream os(file);
	cereal::JSONInputArchive archive(os);
	archive(cereal::make_nvp("springls", *this));
	loadResources(aly::RemoveTrailingSlash(aly::GetParentDirectory(file)));
}
void SpringLevelSet::relax(int iters) {
	Relax<openvdb::util::NullInterrupter> relax(*this);
	for (int iter = 0; iter < iters; iter++) {
		relax.process();
	}
}
void SpringLevelSet::evolve() {
	updateGradient();
	VelocityField grad(*mGradient);
	AdvectionTool advect(*mSignedLevelSet, grad);
	advect.setSpatialScheme(openvdb::math::FIRST_BIAS);
	advect.setTemporalScheme(openvdb::math::TVD_RK2);
	advect.setTrackerSpatialScheme(openvdb::math::FIRST_BIAS);
	advect.setTrackerTemporalScheme(openvdb::math::TVD_RK2);
	advect.advect(0.0, 4.0);
}
template <typename GridType, typename MeshDataAdapter, typename Interrupter>
inline typename GridType::Ptr
meshToVolume(
  Interrupter& interrupter,
  const MeshDataAdapter& mesh,
  SIndexPtr& indexGrid,
  const math::Transform& transform,
  float exteriorBandWidth,
  float interiorBandWidth,
  int flags,
  typename GridType::template ValueConverter<Int32>::Type * polygonIndexGrid=nullptr)
{
    typedef typename GridType::Ptr              GridTypePtr;
    typedef typename GridType::TreeType         TreeType;
    typedef typename TreeType::LeafNodeType     LeafNodeType;
    typedef typename GridType::ValueType        ValueType;

    typedef typename GridType::template ValueConverter<Int32>::Type  Int32GridType;
    typedef typename Int32GridType::TreeType                         Int32TreeType;

    typedef typename TreeType::template ValueConverter<bool>::Type   BoolTreeType;

    //////////

    // Setup

    GridTypePtr distGrid(new GridType(std::numeric_limits<ValueType>::max()));
    distGrid->setTransform(transform.copy());

    ValueType exteriorWidth = ValueType(exteriorBandWidth);
    ValueType interiorWidth = ValueType(interiorBandWidth);

    // Note: inf interior width is all right, this value makes the converter fill
    // interior regions with distance values.
    if (!boost::math::isfinite(exteriorWidth) || boost::math::isnan(interiorWidth)) {
        std::stringstream msg;
        msg << "Illegal narrow band width: exterior = " << exteriorWidth
            << ", interior = " << interiorWidth;
        OPENVDB_LOG_DEBUG(msg.str());
        return distGrid;
    }

    const ValueType voxelSize = ValueType(transform.voxelSize()[0]);

    if (!boost::math::isfinite(voxelSize) || math::isZero(voxelSize)) {
        std::stringstream msg;
        msg << "Illegal transform, voxel size = " << voxelSize;
        OPENVDB_LOG_DEBUG(msg.str());
        return distGrid;
    }

    // Convert narrow band width from voxel units to world space units.
    exteriorWidth *= voxelSize;
    // Avoid the unit conversion if the interior band width is set to
    // inf or std::numeric_limits<float>::max().
    if (interiorWidth < std::numeric_limits<ValueType>::max()) {
        interiorWidth *= voxelSize;
    }

    const bool computeSignedDistanceField = (flags & UNSIGNED_DISTANCE_FIELD) == 0;
    const bool removeIntersectingVoxels = (flags & DISABLE_INTERSECTING_VOXEL_REMOVAL) == 0;
    const bool renormalizeValues = (flags & DISABLE_RENORMALIZATION) == 0;
    const bool trimNarrowBand = (flags & DISABLE_NARROW_BAND_TRIMMING) == 0;
    indexGrid.reset(new Int32GridType(Int32(util::INVALID_IDX)));
    indexGrid->newTree();
    indexGrid->setTransform(transform.copy());

    if (computeSignedDistanceField) {
        distGrid->setGridClass(GRID_LEVEL_SET);
    } else {
        distGrid->setGridClass(GRID_UNKNOWN);
        interiorWidth = ValueType(0.0);
    }

    TreeType& distTree = distGrid->tree();
    Int32TreeType& indexTree = indexGrid->tree();


    // Voxelize mesh

    {
        typedef mesh_to_volume_internal::VoxelizationData<TreeType> VoxelizationDataType;
        typedef tbb::enumerable_thread_specific<typename VoxelizationDataType::Ptr> DataTable;

        DataTable data;
        typedef mesh_to_volume_internal::VoxelizePolygons<TreeType, MeshDataAdapter, Interrupter> Voxelizer;

        const tbb::blocked_range<size_t> polygonRange(0, mesh.polygonCount());

        tbb::parallel_for(polygonRange, Voxelizer(data, mesh, &interrupter));

        for (typename DataTable::iterator i = data.begin(); i != data.end(); ++i) {
            VoxelizationDataType& dataItem = **i;
            mesh_to_volume_internal::combineData(distTree, indexTree, dataItem.distTree, dataItem.indexTree);
        }
    }
    // The progress estimates are based on the observed average time for a few different
    // test cases and is only intended to provide some rough progression feedback to the user.
    if (interrupter.wasInterrupted(30)) return distGrid;


    //////////

    // Classify interior and exterior regions

    if (computeSignedDistanceField) {

        // Determines the inside/outside state for the narrow band of voxels.
        traceExteriorBoundaries(distTree);

        std::vector<LeafNodeType*> nodes;
        nodes.reserve(distTree.leafCount());
        distTree.getNodes(nodes);

        const tbb::blocked_range<size_t> nodeRange(0, nodes.size());

        typedef mesh_to_volume_internal::ComputeIntersectingVoxelSign<TreeType, MeshDataAdapter> SignOp;

        tbb::parallel_for(nodeRange, SignOp(nodes, distTree, indexTree, mesh));

        if (interrupter.wasInterrupted(45)) return distGrid;

        // Remove voxels created by self intersecting portions of the mesh.
        if (removeIntersectingVoxels) {

            tbb::parallel_for(nodeRange,
                mesh_to_volume_internal::ValidateIntersectingVoxels<TreeType>(distTree, nodes));

            tbb::parallel_for(nodeRange,
                mesh_to_volume_internal::RemoveSelfIntersectingSurface<TreeType>(
                    nodes, distTree, indexTree));

            tools::pruneInactive(distTree,  /*threading=*/true);
            tools::pruneInactive(indexTree, /*threading=*/true);
        }
    }

    if (interrupter.wasInterrupted(50)) return distGrid;

    if (distTree.activeVoxelCount() == 0) {
        distGrid.reset((new GridType(ValueType(0.0))));
        return distGrid;
    }

    // Transform values (world space scaling etc.).
    {
        std::vector<LeafNodeType*> nodes;
        nodes.reserve(distTree.leafCount());
        distTree.getNodes(nodes);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, nodes.size()),
            mesh_to_volume_internal::TransformValues<TreeType>(
                nodes, voxelSize, !computeSignedDistanceField));
    }

    // Propagate sign information into tile regions.
    if (computeSignedDistanceField) {
        distTree.root().setBackground(exteriorWidth, /*updateChildNodes=*/false);
        tools::signedFloodFillWithValues(distTree, exteriorWidth, -interiorWidth);
    } else {
        tools::changeBackground(distTree, exteriorWidth);
    }

    if (interrupter.wasInterrupted(54)) return distGrid;


    //////////

    // Expand the narrow band region

    const ValueType minBandWidth = voxelSize * ValueType(2.0);

    if (interiorWidth > minBandWidth || exteriorWidth > minBandWidth) {

        // Create the initial voxel mask.
        BoolTreeType maskTree(false);

        {
            std::vector<LeafNodeType*> nodes;
            nodes.reserve(distTree.leafCount());
            distTree.getNodes(nodes);

            mesh_to_volume_internal::ConstructVoxelMask<TreeType> op(maskTree, distTree, nodes);
            tbb::parallel_reduce(tbb::blocked_range<size_t>(0, nodes.size()), op);
        }

        // Progress estimation
        unsigned maxIterations = std::numeric_limits<unsigned>::max();

        float progress = 54.0f, step = 0.0f;
        double estimated =
            2.0 * std::ceil((std::max(interiorWidth, exteriorWidth) - minBandWidth) / voxelSize);

        if (estimated < double(maxIterations)) {
            maxIterations = unsigned(estimated);
            step = 40.0f / float(maxIterations);
        }

        std::vector<typename BoolTreeType::LeafNodeType*> maskNodes;

        unsigned count = 0;
        while (true) {

            if (interrupter.wasInterrupted(int(progress))) return distGrid;

            const size_t maskNodeCount = maskTree.leafCount();
            if (maskNodeCount == 0) break;

            maskNodes.clear();
            maskNodes.reserve(maskNodeCount);
            maskTree.getNodes(maskNodes);

            const tbb::blocked_range<size_t> range(0, maskNodes.size());

            tbb::parallel_for(range,
                mesh_to_volume_internal::DiffLeafNodeMask<TreeType>(distTree, maskNodes));

            mesh_to_volume_internal::expandNarrowband(distTree, indexTree, maskTree, maskNodes,
                mesh, exteriorWidth, interiorWidth, voxelSize);

            if ((++count) >= maxIterations) break;
            progress += step;
        }
    }

    if (interrupter.wasInterrupted(94)) return distGrid;


    /////////

    // Renormalize distances to smooth out bumps caused by self intersecting
    // and overlapping portions of the mesh and renormalize the level set.

    if (computeSignedDistanceField && renormalizeValues) {

        std::vector<LeafNodeType*> nodes;
        nodes.reserve(distTree.leafCount());
        distTree.getNodes(nodes);

        boost::scoped_array<ValueType> buffer(new ValueType[LeafNodeType::SIZE * nodes.size()]);

        const ValueType offset = ValueType(0.8 * voxelSize);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, nodes.size()),
            mesh_to_volume_internal::OffsetValues<TreeType>(nodes, -offset));

        tbb::parallel_for(tbb::blocked_range<size_t>(0, nodes.size()),
            mesh_to_volume_internal::Renormalize<TreeType>(
                distTree, nodes, buffer.get(), voxelSize));

        tbb::parallel_for(tbb::blocked_range<size_t>(0, nodes.size()),
            mesh_to_volume_internal::MinCombine<TreeType>(nodes, buffer.get()));

        tbb::parallel_for(tbb::blocked_range<size_t>(0, nodes.size()),
            mesh_to_volume_internal::OffsetValues<TreeType>(
                nodes, offset - mesh_to_volume_internal::Tolerance<ValueType>::epsilon()));
    }

    if (interrupter.wasInterrupted(99)) return distGrid;


    /////////

    // Remove active voxels that exceed the narrow band limits

    if (trimNarrowBand && std::min(interiorWidth, exteriorWidth) < voxelSize * ValueType(4.0)) {

        std::vector<LeafNodeType*> nodes;
        nodes.reserve(distTree.leafCount());
        distTree.getNodes(nodes);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, nodes.size()),
            mesh_to_volume_internal::InactivateValues<TreeType>(
                nodes, exteriorWidth, computeSignedDistanceField ? interiorWidth : exteriorWidth));

        tools::pruneLevelSet(
            distTree, exteriorWidth, computeSignedDistanceField ? -interiorWidth : -exteriorWidth);
    }

    return distGrid;
}
template<typename GridType, typename Interrupter>
inline typename boost::enable_if<boost::is_floating_point<typename GridType::ValueType>,
typename GridType::Ptr>::type
doMeshConversion(
    Interrupter& interrupter,
    const openvdb::math::Transform& xform,
    const std::vector<Vec3s>& points,
    const std::vector<Vec3I>& triangles,
    const std::vector<Vec4I>& quads,
	SIndexPtr& indexGrid,
    float exBandWidth,
    float inBandWidth,
    bool unsignedDistanceField = false)
{
    if (points.empty()) {
        return typename GridType::Ptr(new GridType(typename GridType::ValueType(exBandWidth)));
    }

    const size_t numPoints = points.size();
    boost::scoped_array<Vec3s> indexSpacePoints(new Vec3s[numPoints]);

    // transform points to local grid index space
    tbb::parallel_for(tbb::blocked_range<size_t>(0, numPoints),
        mesh_to_volume_internal::TransformPoints<Vec3s>(
                &points[0], indexSpacePoints.get(), xform));

    const int conversionFlags = unsignedDistanceField ? UNSIGNED_DISTANCE_FIELD : 0;

    if (quads.empty()) {

        QuadAndTriangleDataAdapter<Vec3s, Vec3I>
            mesh(indexSpacePoints.get(), numPoints, &triangles[0], triangles.size());

        return ::meshToVolume<GridType>(interrupter,mesh, indexGrid,xform, exBandWidth, inBandWidth, conversionFlags);

    } else if (triangles.empty()) {

        QuadAndTriangleDataAdapter<Vec3s, Vec4I>
            mesh(indexSpacePoints.get(), numPoints, &quads[0], quads.size());

        return ::meshToVolume<GridType>(interrupter,mesh, indexGrid,xform, exBandWidth, inBandWidth, conversionFlags);
    }

    // pack primitives

    const size_t numPrimitives = triangles.size() + quads.size();
    boost::scoped_array<Vec4I> prims(new Vec4I[numPrimitives]);

    for (size_t n = 0, N = triangles.size(); n < N; ++n) {
        const Vec3I& triangle = triangles[n];
        Vec4I& prim = prims[n];
        prim[0] = triangle[0];
        prim[1] = triangle[1];
        prim[2] = triangle[2];
        prim[3] = util::INVALID_IDX;
    }

    const size_t offset = triangles.size();
    for (size_t n = 0, N = quads.size(); n < N; ++n) {
        prims[offset + n] = quads[n];
    }

    QuadAndTriangleDataAdapter<Vec3s, Vec4I>
        mesh(indexSpacePoints.get(), numPoints, prims.get(), numPrimitives);

    return meshToVolume<GridType>(interrupter, mesh, indexGrid,xform, exBandWidth, inBandWidth, conversionFlags);
}
void SpringLevelSet::updateUnSignedLevelSet(double distance) {
	openvdb::math::Transform::Ptr trans = openvdb::math::Transform::createLinearTransform(1.0f);
	using namespace openvdb::tools;
	using namespace openvdb;
    util::NullInterrupter nullInterrupter;
    mUnsignedLevelSet=::doMeshConversion<openvdb::FloatGrid>(nullInterrupter,*trans,mConstellation.mVertexes,mConstellation.mTriangles,  mConstellation.mQuads,mSpringlIndexGrid, distance, distance, true);
	openvdb::tree::LeafManager<FloatGrid::TreeType> lm(mUnsignedLevelSet->tree());
	openvdb::tools::changeBackground(lm, distance);
}


double SpringLevelSet::distanceToConstellation(const Vec3s& pt) {
	openvdb::math::DenseStencil<openvdb::Int32Grid> stencil = openvdb::math::DenseStencil<openvdb::Int32Grid>(*mSpringlIndexGrid, ceil(FILL_DISTANCE));
	stencil.moveTo(Coord((int) floor(pt[0] + 0.5f), (int) floor(pt[1] + 0.5f), (int) floor(pt[2] + 0.5f)));
	int sz = stencil.size();
	double levelSetValue = std::numeric_limits<float>::max();
	std::vector<Index32> stencilCopy;
	stencilCopy.clear();
	Index32 last = -1;
	Index32 springlsCount = mConstellation.getNumSpringls();
	for (int nn = 0; nn < sz; nn++) {
		openvdb::Index32 id = stencil.getValue(nn);
		if (id >= springlsCount)
			continue;
		stencilCopy.push_back(id);
	}
	sz = stencilCopy.size();
	sort(stencilCopy.begin(), stencilCopy.end());
	for (Index32 id : stencilCopy) {
		if (last != id) {
			float d = mConstellation.springls[id].distanceToFaceSqr(pt);
			if (d < levelSetValue) {
				levelSetValue = d;
			}
		}
		last = id;
	}
	return std::sqrt(levelSetValue);
}
void SpringLevelSet::updateSignedLevelSet() {
	openvdb::math::Transform::Ptr trans = openvdb::math::Transform::createLinearTransform(1.0);
	using namespace openvdb::tools;
	using namespace openvdb;
	mSignedLevelSet = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(*trans,mIsoSurface.mVertexes,mIsoSurface.mTriangles, mIsoSurface.mQuads, float(LEVEL_SET_HALF_WIDTH));
}

void SpringLevelSet::updateGradient() {
	//mGradient = openvdb::tools::mGradient(*mUnsignedLevelSet);
	mGradient = advectionForce(*mUnsignedLevelSet);

}
std::list<SpringlNeighbor>& SpringLevelSet::getNearestNeighbors(openvdb::Index32 id, int8_t e) {
	return mNearestNeighbors[mConstellation.springls[id].offset + e];
}
void SpringLevelSet::create(const Constellation& mesh, openvdb::math::Transform::Ptr _transform) {
	this->mTransform = _transform;
	openvdb::math::Transform::Ptr trans = openvdb::math::Transform::createLinearTransform(1.0);
	mSignedLevelSet = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(*trans,mesh.mVertexes,mesh.mTriangles, mesh.mQuads, float(LEVEL_SET_HALF_WIDTH));
	ConvertLevelSetToConstellation(mSignedLevelSet, mIsoSurface);
	mConstellation.copyFrom(mIsoSurface);
	updateIsoSurface();
	for (int iter = 0; iter < 2; iter++) {
		updateUnSignedLevelSet();
		updateNearestNeighbors();
		relax(10);
		updateUnSignedLevelSet(2.5 * openvdb::LEVEL_SET_HALF_WIDTH);
		clean();
		updateUnSignedLevelSet();
		fill();
		fillWithNearestNeighbors();
	}
	updateGradient();

}

int SpringLevelSet::clean() {
	openvdb::math::BoxStencil<openvdb::FloatGrid> stencil(*mSignedLevelSet);
	Vec3s pt, pt1, pt2, pt3;
	float area;
	float minEdgeLength, maxEdgeLength;
	std::vector<Index32> keepList;
	Index32 newVertexCount = 0;
	Index32 newSpringlCount = 0;
	int N = mConstellation.getNumSpringls();
	keepList.reserve(N);
	Index32 index = 0;
	double minls = 1E30, bias = 0, maxls = -1E30, meanls = 0;
	int count = 0;
	int removeFarCount = 0;
	int removeSmallCount = 0;
	int removeAspectCount = 0;
	std::vector<float> levelSetValues(mConstellation.springls.size());
#pragma omp for
	for (int nn = 0; nn < (int) mConstellation.springls.size(); nn++) {
		Springl& springl = mConstellation.springls[nn];
		pt = springl.particle();
		stencil.moveTo(Coord(std::floor(pt[0]), std::floor(pt[1]), std::floor(pt[2])));
		levelSetValues[nn] = stencil.interpolation(pt);
	}
	for (Springl& springl : mConstellation.springls) {
		float levelSetValue = levelSetValues[count++];
		int K = springl.size();
		double v = std::abs(levelSetValue);
		meanls += v;
		bias += levelSetValue;
		minls = std::min(minls, v);
		maxls = std::max(maxls, v);
		if (v <= CLEAN_DISTANCE) {
			minEdgeLength = 1E30;
			maxEdgeLength = -1E30;
			area = 0.0f;
			for (int i = 0; i < K; i++) {
				pt1 = springl[i];
				pt2 = springl[(i + 1) % K];
				float len = (pt1 - pt2).length();
				minEdgeLength = std::min(minEdgeLength, len);
				maxEdgeLength = std::max(maxEdgeLength, len);
			}
			float aspect = minEdgeLength / maxEdgeLength;
			area = springl.area();
			//std::cout<<"AREA "<<area<<" ASPECT "<<aspect<<std::endl;
			if (area >= MIN_AREA && area < MAX_AREA && aspect >= MIN_ASPECT_RATIO) {
				keepList.push_back(springl.id);
				newSpringlCount++;
				newVertexCount += K;
			} else {
				if (area < MIN_AREA || area >= MAX_AREA) {
					removeSmallCount++;
				}
				if (aspect < MIN_ASPECT_RATIO) {
					removeAspectCount++;
				}
			}
		} else {
			removeFarCount++;
		}
		index++;
	}
	levelSetValues.clear();
	meanls /= count;
	bias /= count;
	std::cout << "Clean: mean=" << meanls << " bias=" << bias << " Level Set= [" << minls << "," << maxls << "] [far:" << removeFarCount << ", small:"<< removeSmallCount << ", aspect:" << removeAspectCount << "]" << std::endl;
	if ((int) newSpringlCount == N)return 0;
	Index32 springlOffset = 0;
	Index32 vertexOffset = 0;
	Index32 quadIndex = 0;
	Index32 triIndex = 0;
	for (int n : keepList) {
		Springl& rspringl = mConstellation.springls[n];
		Springl& springl = mConstellation.springls[springlOffset];
		int K = rspringl.size();
		if ((int) springlOffset != n) {
			mConstellation.mParticles[springlOffset] = mConstellation.mParticles[n];

			if (mConstellation.mParticleVelocity.size() > 0) {
				mConstellation.mParticleVelocity[springlOffset] = mConstellation.mParticleVelocity[n];
			}
			if (mConstellation.mVertexVelocity.size() > 0) {
				for (int n = 0; n < springl.size(); n++) {
					mConstellation.mVertexVelocity[springl.offset + n] = mConstellation.mVertexVelocity[rspringl.offset + n];
				}
			}
			if (mConstellation.mParticleLabel.size() > 0) {
				mConstellation.mParticleLabel[springlOffset] = mConstellation.mParticleLabel[n];
			}
			mConstellation.mParticleNormals[springlOffset] = mConstellation.mParticleNormals[n];
			springl.offset = vertexOffset;
			springl.id = springlOffset;

			for (int k = 0; k < K; k++) {
				mConstellation.mVertexes[vertexOffset + k] = mConstellation.mVertexes[rspringl.offset + k];
				mConstellation.mVertexNormals[vertexOffset + k] = mConstellation.mVertexNormals[rspringl.offset + k];
			}
			if (K == 4) {
				Vec4I quad;
				for (int k = 0; k < K; k++) {
					quad[k] = vertexOffset + k;
				}
				mConstellation.mQuads[quadIndex++]=quad;
			} else if (K == 3) {
				Vec3I tri;
				for (int k = 0; k < K; k++) {
					tri[k] = vertexOffset + k;
				}
				mConstellation.mTriangles[triIndex++]=tri;
			}

		} else {
			if (K == 4) {
				quadIndex ++;
			} else if (K == 3) {
				triIndex ++;
			}
		}
		vertexOffset += K;
		springlOffset++;
	}
	mConstellation.mQuads.erase(mConstellation.mQuads.begin() + quadIndex, mConstellation.mQuads.end());
	mConstellation.mTriangles.erase(mConstellation.mTriangles.begin() + triIndex, mConstellation.mTriangles.end());
	mConstellation.springls.erase(mConstellation.springls.begin() + springlOffset, mConstellation.springls.end());
	mConstellation.mParticles.erase(mConstellation.mParticles.begin() + springlOffset, mConstellation.mParticles.end());

	if (mConstellation.mParticleVelocity.size() > 0) {
		mConstellation.mParticleVelocity.erase(mConstellation.mParticleVelocity.begin() + springlOffset, mConstellation.mParticleVelocity.end());
	}
	if (mConstellation.mVertexVelocity.size() > 0) {
		mConstellation.mVertexVelocity.erase(mConstellation.mVertexVelocity.begin() + vertexOffset, mConstellation.mVertexVelocity.end());
	}
	mConstellation.mParticleNormals.erase(mConstellation.mParticleNormals.begin() + springlOffset, mConstellation.mParticleNormals.end());
	mConstellation.mVertexNormals.erase(mConstellation.mVertexNormals.begin() + vertexOffset, mConstellation.mVertexNormals.end());
	mConstellation.mVertexes.erase(mConstellation.mVertexes.begin() + vertexOffset, mConstellation.mVertexes.end());
	mCleanCount += (N - newSpringlCount);
	return (N - newSpringlCount);
}
void SpringLevelSet::create(FloatGrid& grid) {
	this->mTransform = grid.transformPtr();
	grid.setTransform(openvdb::math::Transform::createLinearTransform(1.0));
	mSignedLevelSet = boost::static_pointer_cast<FloatGrid>(grid.copyGrid(CopyPolicy::CP_COPY));
	ConvertLevelSetToConstellation(mSignedLevelSet, mIsoSurface);
	updateSignedLevelSet();
	mConstellation.copyFrom(mIsoSurface);
	updateIsoSurface();
	for (int iter = 0; iter < 2; iter++) {
		updateUnSignedLevelSet();
		updateNearestNeighbors();
		relax(10);
		updateUnSignedLevelSet(2.5 * openvdb::LEVEL_SET_HALF_WIDTH);
		clean();
		updateUnSignedLevelSet();
		fill();
		fillWithNearestNeighbors();
	}
	updateGradient();
}
void SpringLevelSet::create(RegularGrid<float>& grid) {
	this->mTransform = grid.transformPtr();
	mSignedLevelSet = std::unique_ptr<FloatGrid>(new FloatGrid());
	openvdb::tree::LeafManager<FloatGrid::TreeType> lm(mSignedLevelSet->tree());
	openvdb::tools::changeBackground(lm, openvdb::LEVEL_SET_HALF_WIDTH);
	mSignedLevelSet->setTransform(grid.transformPtr());
	openvdb::tools::copyFromDense(grid, *mSignedLevelSet, 0.25f);
	mSignedLevelSet->setTransform(Transform::createLinearTransform(1.0));
	ConvertLevelSetToConstellation(mSignedLevelSet, mIsoSurface);
	updateSignedLevelSet();
	mConstellation.copyFrom(mIsoSurface);
	updateIsoSurface();
	for (int iter = 0; iter < 2; iter++) {
		updateUnSignedLevelSet();
		updateNearestNeighbors();
		relax(10);
		updateUnSignedLevelSet(2.5 * openvdb::LEVEL_SET_HALF_WIDTH);
		clean();
		updateUnSignedLevelSet();
		fill();
		fillWithNearestNeighbors();
	}
	updateGradient();
}
void SpringLevelSet::saveResources(const std::string& directory) {
	constellationFile = aly::MakeString() << directory << ALY_PATH_SEPARATOR<<"constellation_"<<std::setw(5)<<std::setfill('0')<<mIteration<<".ply";
	isoSurfaceFile=aly::MakeString()<<directory<<ALY_PATH_SEPARATOR<<"iso-surface_"<<std::setw(5)<<std::setfill('0')<<mIteration<<".ply";
	signedLevelSetFile=aly::MakeString()<<directory<<ALY_PATH_SEPARATOR<<"signed_levelset_"<<std::setw(5)<<std::setfill('0')<<mIteration<<".vdb";
	WriteConstellationToFile(constellationFile,mConstellation);
	WriteConstellationToFile(isoSurfaceFile,mIsoSurface);
	/*
	try {
		openvdb::io::File file(signedLevelSetFile);
		openvdb::GridPtrVec grids;
		mSignedLevelSet->transform()=transform();
		grids.push_back(mSignedLevelSet);
		file.write(grids);
		mSignedLevelSet->setTransform(Transform::createLinearTransform(1.0));
	} catch(openvdb::Exception& e) {
		std::cerr<<"OpenVDB Error: "<<e.what()<<std::endl;
	} catch(std::exception& e) {
		std::cerr<<e.what()<<std::endl;
	}
	*/
}
void SpringLevelSet::loadResources(const std::string& directory) {
	ReadConstellationFromFile(constellationFile,mConstellation);
	ReadConstellationFromFile(isoSurfaceFile,mIsoSurface);
	/*
	try {
		openvdb::io::File file(signedLevelSetFile);
		file.open();
		if(file.isOpen()){
			openvdb::GridPtrVecPtr grids =file.getGrids();
			openvdb::GridPtrVec allGrids;
			allGrids.insert(allGrids.end(), grids->begin(), grids->end());
			GridBase::Ptr ptr = allGrids[0];
			mSignedLevelSet=boost::static_pointer_cast<FloatGrid>(ptr);
			transform()=mSignedLevelSet->transform();
			mSignedLevelSet->setTransform(openvdb::math::Transform::createLinearTransform(1.0));
		}
	} catch(openvdb::Exception& e) {
		std::cerr<<"OpenVDB Error: "<<e.what()<<std::endl;
	} catch(std::exception& e) {
		std::cerr<<e.what()<<std::endl;
	}
	*/
}
void SpringLevelSet::updateIsoSurface() {
	mVolToMesh(*mSignedLevelSet);
	ConvertLevelSetToConstellation(mSignedLevelSet, mVolToMesh, mIsoSurface);
}
int SpringLevelSet::fill() {

	openvdb::tree::ValueAccessor<FloatGrid::TreeType> acc(mSignedLevelSet->tree());
	openvdb::math::GenericMap map(mSignedLevelSet->transform());

	openvdb::math::DenseStencil<openvdb::Int32Grid> stencil = openvdb::math::DenseStencil<openvdb::Int32Grid>(*mSpringlIndexGrid, ceil(FILL_DISTANCE));

	openvdb::tools::PolygonPoolList& polygonPoolList = mVolToMesh.polygonPoolList();
	Vec3s p[4];
	Vec3s refPoint;

	Index32 springlsCount = mConstellation.getNumSpringls();
	Index32 pcounter = mConstellation.getNumSpringls();
	Index32 counter = mConstellation.getNumVertexes();
	int added = 0;

	const float D2 = FILL_DISTANCE * FILL_DISTANCE;
	Index64 N = mVolToMesh.polygonPoolListSize();
	Index64 I;
	fillList.clear();
	for (Index64 n = 0; n < N; ++n) {
		const openvdb::tools::PolygonPool& polygons = polygonPoolList[n];
		I = polygons.numQuads();
#pragma omp for
		for (Index64 i = 0; i < I; ++i) {
			std::vector<openvdb::Index32> stencilCopy;
			const openvdb::Vec4I& quad = polygons.quad(i);
			PointList& pl=mVolToMesh.pointList();
			p[0] = pl[quad[3]];
			p[1] = pl[quad[2]];
			p[2] = pl[quad[1]];
			p[3] = pl[quad[0]];
			refPoint = 0.25f * (p[0] + p[1] + p[2] + p[3]);
			stencil.moveTo(Coord(std::floor(refPoint[0] + 0.5f), std::floor(refPoint[1] + 0.5f), std::floor(refPoint[2] + 0.5f)));
			int sz = stencil.size();
			float levelSetValue = std::numeric_limits<float>::max();
			int last = -1;
			for (int nn = 0; nn < sz; nn++) {
				openvdb::Index32 id = stencil.getValue(nn);
				if (id >= springlsCount)
					continue;
				stencilCopy.push_back(id);
			}
			sz = stencilCopy.size();
			sort(stencilCopy.begin(), stencilCopy.end());
			for (int nn = 0; nn < sz; nn++) {
				openvdb::Index32 id = stencilCopy[nn];
				if (last != (int) id) {
					float d = mConstellation.springls[id].distanceToFaceSqr(refPoint);
					if (d < levelSetValue) {
						levelSetValue = d;
					}
				}
				last = id;
			}
			if (levelSetValue > D2) {
#pragma omp critical
				{
					added++;
					mConstellation.mVertexes.push_back(p[0]);
					mConstellation.mVertexes.push_back(p[1]);
					mConstellation.mVertexes.push_back(p[2]);
					mConstellation.mVertexes.push_back(p[3]);
					Springl springl(&mConstellation);
					springl.offset = counter;
					springl.id = mConstellation.springls.size();
					mConstellation.mQuads.push_back(Vec4I(counter, counter + 1, counter + 2, counter + 3));
					mConstellation.mParticles.push_back(springl.computeCentroid());
					if (mConstellation.mParticleVelocity.size() > 0) {
						fillList.push_back(springl.id);
						mConstellation.mParticleVelocity.push_back(Vec3s(0.0f));
					}
					if (mConstellation.mVertexVelocity.size() > 0) {
						mConstellation.mVertexVelocity.push_back(Vec3s(0.0f));
						mConstellation.mVertexVelocity.push_back(Vec3s(0.0f));
						mConstellation.mVertexVelocity.push_back(Vec3s(0.0f));
						mConstellation.mVertexVelocity.push_back(Vec3s(0.0f));
					}
					if (mConstellation.mParticleLabel.size() > 0) {
						mConstellation.mParticleLabel.push_back(0);
					}
					openvdb::Vec3s norm = springl.computeNormal();
					mConstellation.mParticleNormals.push_back(norm);
					mConstellation.mVertexNormals.push_back(norm);
					mConstellation.mVertexNormals.push_back(norm);
					mConstellation.mVertexNormals.push_back(norm);
					mConstellation.mVertexNormals.push_back(norm);
					mConstellation.springls.push_back(springl);
					pcounter++;
					counter += 4;
				}
			}
		}
		I = polygons.numTriangles();
#pragma omp for
		for (Index64 i = 0; i < I; ++i) {
			std::vector<openvdb::Index32> stencilCopy;
			const openvdb::Vec3I& tri = polygons.triangle(i);
			PointList& pl=mVolToMesh.pointList();
			p[0] = pl[tri[2]];
			p[1] = pl[tri[1]];
			p[2] = pl[tri[0]];
			refPoint = 0.25f * (p[0] + p[1] + p[2]);
			stencil.moveTo(Coord(std::floor(refPoint[0] + 0.5f), std::floor(refPoint[1] + 0.5f), std::floor(refPoint[2] + 0.5f)));
			int sz = stencil.size();
			float levelSetValue = std::numeric_limits<float>::max();
			int last = -1;
			for (int nn = 0; nn < sz; nn++) {
				openvdb::Index32 id = stencil.getValue(nn);
				if (id >= springlsCount)
					continue;
				stencilCopy.push_back(id);
			}
			sz = stencilCopy.size();
			sort(stencilCopy.begin(), stencilCopy.end());
			for (int nn = 0; nn < sz; nn++) {
				openvdb::Index32 id = stencilCopy[nn];
				if (last != (int) id) {
					float d = mConstellation.springls[id].distanceToFaceSqr(refPoint);
					if (d < levelSetValue) {
						levelSetValue = d;
					}
				}
				last = id;
			}
			if (levelSetValue > D2) {
#pragma omp critical
				{
					added++;
					mConstellation.mVertexes.push_back(p[0]);
					mConstellation.mVertexes.push_back(p[1]);
					mConstellation.mVertexes.push_back(p[2]);
					Springl springl(&mConstellation);
					springl.offset = counter;
					springl.id = mConstellation.springls.size();

					mConstellation.mTriangles.push_back(Vec3I(counter, counter + 1, counter + 2));
					mConstellation.mParticles.push_back(springl.computeCentroid());
					if (mConstellation.mParticleVelocity.size() > 0) {
						fillList.push_back(springl.id);
						mConstellation.mParticleVelocity.push_back(Vec3s(0.0f));
					}
					if (mConstellation.mVertexVelocity.size() > 0) {
						mConstellation.mVertexVelocity.push_back(Vec3s(0.0f));
						mConstellation.mVertexVelocity.push_back(Vec3s(0.0f));
						mConstellation.mVertexVelocity.push_back(Vec3s(0.0f));
					}
					if (mConstellation.mParticleLabel.size() > 0) {
						mConstellation.mParticleLabel.push_back(0);
					}
					openvdb::Vec3s norm = springl.computeNormal();

					mConstellation.mParticleNormals.push_back(norm);
					mConstellation.mVertexNormals.push_back(norm);
					mConstellation.mVertexNormals.push_back(norm);
					mConstellation.mVertexNormals.push_back(norm);
					mConstellation.springls.push_back(springl);
					pcounter++;
					counter += 3;
				}
			}
		}
	}

	mFillCount += added;
	return added;
}
void SpringLevelSet::fillWithVelocityField(MACGrid<float>& grid, float radius) {
	for (int fid : fillList) {
		Springl& springl = mConstellation.springls[fid];
		//Is particle() in the correct coordinate space?
		Vec3d wpt = transform().indexToWorld(springl.particle());
		mConstellation.mParticleVelocity[fid] = grid.maxInterpolate(Vec3s(wpt), radius);
		for (int n = 0; n < springl.size(); n++) {
			wpt = transform().indexToWorld(springl[n]);
			mConstellation.mVertexVelocity[springl.offset + n] = grid.maxInterpolate(Vec3s(wpt), radius);
		}
	}
	fillList.clear();
}

void SpringLevelSet::fillWithNearestNeighbors() {
	if (fillList.size() > 0) {
		updateUnSignedLevelSet();
		updateNearestNeighbors();
		for (int cycle = 0; cycle < 16; cycle++) {
			int unfilledCount = 0;
			for (int fid : fillList) {
				Springl& springl = mConstellation.springls[fid];
				int K = springl.size();
				Vec3s vel(0.0);
				double wsum = 0.0;
				for (int k = 0; k < K; k++) {
					std::list<SpringlNeighbor>& map = getNearestNeighbors(springl.id, k);
					for (SpringlNeighbor ci : map) {
						Springl& nbr = getSpringl(ci.springlId);
						Vec3s v = nbr.particleVelocity();
						if (v.lengthSqr() > 0) {
							vel += v;
							wsum += 1.0f;
						}
					}
				}
				if (wsum > 0.0f) {
					vel *= 1.0 / wsum;
					mConstellation.mParticleVelocity[fid] = vel;
					for (int n = 0; n < springl.size(); n++) {
						mConstellation.mVertexVelocity[springl.offset + n] = vel;
					}
				} else {
					unfilledCount++;
				}
			}
			if (unfilledCount == 0)
				break;
			std::cout << cycle << ":: un-filled " << unfilledCount << std::endl;
		}
		fillList.clear();
	}
}
void SpringLevelSet::computeStatistics(const Constellation& mesh) {
	std::vector<Index32> keepList;
	int N = mConstellation.getNumSpringls();
	keepList.reserve(N);
	double minls = 1E30, bias = 0, maxls = -1E30, meanls = 0, v, sqrs = 0, stdev;
	int count = 0;
	for (Vec3s pt : mesh.mVertexes) {
		float levelSetValue = distanceToConstellation(pt);
		count++;
		v = fabs(levelSetValue);
		sqrs += v * v;
		meanls += v;
		bias += levelSetValue;
		minls = std::min(minls, v);
		maxls = std::max(maxls, v);
		count++;
	}
	meanls /= count;
	bias /= count;
	stdev = std::sqrt(sqrs / count - meanls * meanls);
	std::cout << ">>Constellation Vertex mean=" << meanls << " std dev.=" << stdev << " bias=" << bias << " [" << minls << "," << maxls << "]" << std::endl;
	meanls = 0;
	bias = 0;
	sqrs = 0;
	minls = 1E30;
	maxls = -1E30;
	count = 0;

	for (Vec3s pt : mesh.mParticles) {
		float levelSetValue = distanceToConstellation(pt);
		count++;
		v = fabs(levelSetValue);
		sqrs += v * v;
		meanls += v;
		bias += levelSetValue;
		minls = std::min(minls, v);
		maxls = std::max(maxls, v);
		count++;
	}
	meanls /= count;
	bias /= count;

	stdev = std::sqrt(sqrs / count - meanls * meanls);
	if (count > 0)
		std::cout << ">>Constellation Particle mean=" << meanls << " std dev.=" << stdev << " bias=" << bias << " [" << minls << "," << maxls << "]"
				<< std::endl;
}
void SpringLevelSet::computeStatistics(const Constellation& mesh, FloatGrid& levelSet) {
	openvdb::math::BoxStencil<openvdb::FloatGrid> stencil(levelSet);
	std::vector<Index32> keepList;
	int N = mConstellation.getNumSpringls();
	keepList.reserve(N);
	double minls = 1E30, bias = 0, maxls = -1E30, meanls = 0, v, sqrs = 0, stdev;
	int count = 0;
	for (Vec3s pt : mesh.mVertexes) {
		stencil.moveTo(Coord(std::floor(pt[0]), std::floor(pt[1]), std::floor(pt[2])));
		float levelSetValue = stencil.interpolation(pt);
		count++;
		v = fabs(levelSetValue);
		sqrs += v * v;
		meanls += v;
		bias += levelSetValue;
		minls = std::min(minls, v);
		maxls = std::max(maxls, v);
		count++;
	}
	meanls /= count;
	bias /= count;

	stdev = std::sqrt(sqrs / count - meanls * meanls);
	std::cout << ">>Vertex mean=" << meanls << " std dev.=" << stdev << " bias=" << bias << " [" << minls << "," << maxls << "]" << std::endl;

	meanls = 0;
	bias = 0;
	sqrs = 0;
	minls = 1E30;
	maxls = -1E30;
	count = 0;

	for (Vec3s pt : mesh.mParticles) {
		stencil.moveTo(Coord(std::floor(pt[0]), std::floor(pt[1]), std::floor(pt[2])));
		float levelSetValue = stencil.interpolation(pt);
		count++;
		v = fabs(levelSetValue);
		sqrs += v * v;
		meanls += v;
		bias += levelSetValue;
		minls = std::min(minls, v);
		maxls = std::max(maxls, v);
		count++;
	}
	meanls /= count;
	bias /= count;

	stdev = std::sqrt(sqrs / count - meanls * meanls);
	if (count > 0)
		std::cout << ">>Particle mean=" << meanls << " std dev.=" << stdev << " bias=" << bias << " [" << minls << "," << maxls << "]" << std::endl;

}
