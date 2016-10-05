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

#ifndef VDBTOOLS_INCLUDE_ADVECTOPERATION_H_
#define VDBTOOLS_INCLUDE_ADVECTOPERATION_H_
#include "SpringLevelSet.h"
#include "MeshVertexRange.h"
#include <openvdb/openvdb.h>
template<typename FieldT> openvdb::Vec3d ComputeVelocity(const FieldT& field,
		TemporalIntegrationScheme scheme, openvdb::Vec3d pt, double t,
		double h);

template<typename FieldT> class AdvectSpringlOperation {
private:
	TemporalIntegrationScheme mIntegrationScheme;
public:
	AdvectSpringlOperation(
			TemporalIntegrationScheme integrationScheme =
					TemporalIntegrationScheme::UNKNOWN_TIS) :
			mIntegrationScheme(integrationScheme) {
	}
	void compute(Springl& springl, SpringLevelSet& mGrid, const FieldT& field,
			double t, double h) {
		using namespace openvdb;
		openvdb::math::Transform::Ptr trans = mGrid.transformPtr();
		Vec3d v = Vec3d(springl.particle());
		Vec3d pt = trans->indexToWorld(v);
		Vec3d vel = ComputeVelocity(field, mIntegrationScheme, pt, t, h);
		springl.particle() = trans->worldToIndex(pt + vel);	//Apply integration scheme here, need buffer for previous time points?
		int K = springl.size();
		for (int k = 0; k < K; k++) {
			pt = trans->indexToWorld(springl[k]);
			vel = ComputeVelocity(field, mIntegrationScheme, pt, t, h);
			springl[k] = trans->worldToIndex(pt + vel);
		}
	}
	double findTimeStep(Springl& springl, SpringLevelSet& mGrid,
			const FieldT& field, double t) {
		using namespace openvdb;
		openvdb::math::Transform::Ptr trans = mGrid.transformPtr();
		Vec3d v = Vec3d(springl.particle());
		Vec3d pt = trans->indexToWorld(v);
		Vec3d vec = field(pt, t);
		return std::max(std::max(fabs(vec[0]), fabs(vec[1])), fabs(vec[2]));
	}

};
template<typename FieldT> class AdvectParticleOperation {
private:
	TemporalIntegrationScheme mIntegrationScheme;
public:
	AdvectParticleOperation(
			TemporalIntegrationScheme integrationScheme =
					TemporalIntegrationScheme::UNKNOWN_TIS) :
			mIntegrationScheme(integrationScheme) {
	}
	void compute(Springl& springl, SpringLevelSet& mGrid, const FieldT& field,
			double t, double h) {
		using namespace openvdb;
		openvdb::math::Transform::Ptr trans = mGrid.transformPtr();
		Vec3d v = Vec3d(springl.particle());
		Vec3d pt = trans->indexToWorld(v);
		Vec3d vel = ComputeVelocity(field, mIntegrationScheme, pt, t, h);
		springl.particle() = trans->worldToIndex(pt + vel);	//Apply integration scheme here, need buffer for previous time points?
	}
	double findTimeStep(Springl& springl, SpringLevelSet& mGrid,
			const FieldT& field, double t) {
		using namespace openvdb;
		openvdb::math::Transform::Ptr trans = mGrid.transformPtr();
		Vec3d v = Vec3d(springl.particle());
		Vec3d pt = trans->indexToWorld(v);
		Vec3d vec = field(pt, t);
		return std::max(std::max(fabs(vec[0]), fabs(vec[1])), fabs(vec[2]));
	}

};
template<typename FieldT> class AdvectParticleAndVertexOperation {
private:
	TemporalIntegrationScheme mIntegrationScheme;
public:
	AdvectParticleAndVertexOperation(
			TemporalIntegrationScheme integrationScheme =
					TemporalIntegrationScheme::UNKNOWN_TIS) :
			mIntegrationScheme(integrationScheme) {
	}
	void compute(Springl& springl, SpringLevelSet& mGrid, const FieldT& field,
			double t, double h) {
		using namespace openvdb;
		openvdb::math::Transform::Ptr trans = mGrid.transformPtr();
		Vec3d v = Vec3d(springl.particle());
		Vec3d pt = trans->indexToWorld(v);
		Vec3d vel = ComputeVelocity(field, mIntegrationScheme, pt, t, h);
		springl.particle() = trans->worldToIndex(pt + vel);	//Apply integration scheme here, need buffer for previous time points?

		for(int k=0;k<springl.size();k++){
			v = Vec3d(springl[k]);
			pt = trans->indexToWorld(v);
			vel = ComputeVelocity(field, mIntegrationScheme, pt, t, h);
			springl[k] = trans->worldToIndex(pt + vel);	//Apply integration scheme here, need buffer for previous time points?
		}
	}
	double findTimeStep(Springl& springl, SpringLevelSet& mGrid,
			const FieldT& field, double t) {
		using namespace openvdb;
		openvdb::math::Transform::Ptr trans = mGrid.transformPtr();
		Vec3d v = Vec3d(springl.particle());
		Vec3d pt = trans->indexToWorld(v);
		Vec3d vec = field(pt, t);
		return std::max(std::max(fabs(vec[0]), fabs(vec[1])), fabs(vec[2]));
	}

};

template<typename FieldT> class AdvectMeshVertexOperation {
private:
	TemporalIntegrationScheme mIntegrationScheme;
public:
	AdvectMeshVertexOperation(
			TemporalIntegrationScheme integrationScheme =
					TemporalIntegrationScheme::UNKNOWN_TIS) :
			mIntegrationScheme(integrationScheme) {
	}
	double findTimeStep(size_t vid, SpringLevelSet& mGrid, Constellation& mMesh,
			const FieldT& field, double t) {
		using namespace openvdb;
		openvdb::math::Transform::Ptr trans = mGrid.transformPtr();
		Vec3s vert = mGrid.mIsoSurface.mVertexes[vid];
		Vec3d v = Vec3d(vert);
		Vec3d pt = trans->indexToWorld(v);
		Vec3f vec = field(pt, t);
		return std::max(std::max(fabs(vec[0]), fabs(vec[1])), fabs(vec[2]));
	}
	void compute(size_t vid, SpringLevelSet& mGrid, const FieldT& field,
			double t, double dt) {
		using namespace openvdb;
		openvdb::math::Transform::Ptr trans = mGrid.transformPtr();
		Vec3s vert = mGrid.mIsoSurface.mVertexes[vid];
		Vec3d v = Vec3d(vert);
		Vec3d pt = trans->indexToWorld(v);
		Vec3d vel = ComputeVelocity(field, mIntegrationScheme, pt, t, dt);
		mGrid.mIsoSurface.mVertexes[vid] = trans->worldToIndex(pt + vel);
	}
};


template<typename OperatorT, typename FieldT,
		typename InterruptT = openvdb::util::NullInterrupter>
class AdvectSpringlFieldOperator {
public:
	SpringLevelSet& mGrid;
	AdvectSpringlFieldOperator(SpringLevelSet& grid, const FieldT& field,
			TemporalIntegrationScheme scheme, double t, double dt,
			InterruptT* _interrupt) :
			mGrid(grid), mField(field), mIntegrationScheme(scheme), mInterrupt(
					_interrupt), mTime(t), mTimeStep(dt) {

	}
	virtual ~AdvectSpringlFieldOperator() {
	}
	void process(bool threaded = true) {
		if (mInterrupt)
			mInterrupt->start("Processing springls");
		SpringlRange range(mGrid.mConstellation);
		if (threaded) {
			tbb::parallel_for(range, *this);
		} else {
			(*this)(range);
		}

		if (mInterrupt)
			mInterrupt->end();
	}

	/// @note Never call this public method directly - it is called by
	/// TBB threads only!
	void operator()(const SpringlRange& range) const {
		if (openvdb::util::wasInterrupted(mInterrupt))
			tbb::task::self().cancel_group_execution();
		OperatorT OpT(mIntegrationScheme);
		for (typename SpringlRange::Iterator springl = range.begin(); springl;
				++springl) {
			OpT.compute(*springl, mGrid, mField, mTime, mTimeStep);
		}
	}

protected:
	double mTime;
	double mTimeStep;
	TemporalIntegrationScheme mIntegrationScheme;
	const FieldT& mField;
	InterruptT* mInterrupt;

};
template<typename OperatorT, typename FieldT,
		typename InterruptT = openvdb::util::NullInterrupter>
class AdvectMeshVertexOperator {
public:
	SpringLevelSet& mGrid;
	AdvectMeshVertexOperator(SpringLevelSet& grid, const FieldT& field,
			TemporalIntegrationScheme scheme, double t, double dt,
			InterruptT* _interrupt) :
			mGrid(grid), mField(field), mInterrupt(_interrupt), mIntegrationScheme(
					scheme), mTime(t), mTimeStep(dt) {

	}
	virtual ~AdvectMeshVertexOperator() {
	}
	void process(bool threaded = true) {
		if (mInterrupt)
			mInterrupt->start("Processing springls");
		MeshVertexRange range(mGrid.mIsoSurface);
		if (threaded) {
			tbb::parallel_for(range, *this);
		} else {
			(*this)(range);
		}

		if (mInterrupt)
			mInterrupt->end();
	}
	/// @note Never call this public method directly - it is called by
	/// TBB threads only!
	void operator()(const MeshVertexRange& range) const {
		if (openvdb::util::wasInterrupted(mInterrupt))
			tbb::task::self().cancel_group_execution();
		OperatorT OpT(mIntegrationScheme);
		for (typename MeshVertexRange::Iterator vert = range.begin(); vert;
				++vert) {
			OpT.compute(*vert, mGrid, mField, mTime, mTimeStep);
		}
	}

protected:
	double mTime;
	double mTimeStep;
	TemporalIntegrationScheme mIntegrationScheme;

	const FieldT& mField;
	InterruptT* mInterrupt;

};
template<typename FieldT> openvdb::Vec3d ComputeVelocity(const FieldT& field,
		TemporalIntegrationScheme scheme, openvdb::Vec3d pt, double t,
		double h) {
	openvdb::Vec3d velocity(0.0);
	openvdb::Vec3d k1, k2, k3, k4;
	switch (scheme) {
	case TemporalIntegrationScheme::RK1:
		velocity = h * field(pt, t);
		break;
	case TemporalIntegrationScheme::RK2:
		k1 = h * field(pt, t);
		velocity = h * field(pt + 0.5 * k1, t + 0.5f * h);
		break;
	case TemporalIntegrationScheme::RK3:
		k1 = h * field(pt, t);
		k2 = h * field(pt + 0.5 * k1, t + 0.5f * h);
		k3 = h * field(pt - 1.0 * k1 + 2.0 * k2, t + h);
		velocity = (1.0f / 6.0f) * (k1 + 4 * k2 + k3);
		break;
	case TemporalIntegrationScheme::RK4a:
		k1 = h * field(pt, t);
		k2 = h * field(pt + 0.5f * k1, t + 0.5f * h);
		k3 = h * field(pt + 0.5f * k2, t + 0.5f * h);
		k4 = h * field(pt + k3, t + h);
		velocity = (1.0f / 6.0f) * (k1 + 2 * k2 + 2 * k3 + k4);
		break;
	case TemporalIntegrationScheme::RK4b:
	default:
		k1 = h * field(pt, t);
		k2 = h * field(pt + (1 / 3.0) * k1, t + (1 / 3.0) * h);
		k3 = h * field(pt - (1 / 3.0) * k1 + k2, t + (2 / 3.0) * h);
		k4 = h * field(pt + k1 - k2 + k3, t + h);
		velocity = (1.0f / 8.0f) * (k1 + 3 * k2 + 3 * k3 + k4);
		break;

	}
	return velocity;
}
#endif /* VDBTOOLS_INCLUDE_ADVECTOPERATION_H_ */
