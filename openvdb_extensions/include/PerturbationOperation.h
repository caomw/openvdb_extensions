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

#ifndef VDBTOOLS_INCLUDE_PERTURBATIONOPERATION_H_
#define VDBTOOLS_INCLUDE_PERTURBATIONOPERATION_H_

#include "SpringLevelSet.h"
#include "SpringlRange.h"
#include "MeshVertexRange.h"
#include <openvdb/util/NullInterrupter.h>
template<typename OperatorT,
		typename InterruptT = openvdb::util::NullInterrupter>
class PerturbSpringlOperator {
protected:
	InterruptT* mInterrupt;

public:
	SpringLevelSet& mGrid;
	double mDt;
	PerturbSpringlOperator(SpringLevelSet& grid, InterruptT* _interrupt,
			double dt) :
				 mInterrupt(_interrupt), mGrid(grid),mDt(dt) {

	}
	virtual ~PerturbSpringlOperator() {
	}
	void process(bool threaded = true) {
		if (mInterrupt)
			mInterrupt->start("Processing springls");
		OperatorT::init(mGrid);
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
		for (typename SpringlRange::Iterator springl = range.begin(); springl;++springl) {
			OperatorT::apply(*springl, mGrid, mDt);
		}
	}


};
template<typename OperatorT,
		typename InterruptT = openvdb::util::NullInterrupter>
class PerturbMeshVertexOperator {
public:
	SpringLevelSet& mGrid;
	double mDt;
	PerturbMeshVertexOperator(SpringLevelSet& grid, InterruptT* _interrupt,
			double dt) :
			mGrid(grid), mInterrupt(_interrupt), mDt(dt) {

	}
	virtual ~PerturbMeshVertexOperator() {
	}
	void process(bool threaded = true) {
		if (mInterrupt)
			mInterrupt->start("Processing springls");
		OperatorT::init(mGrid);
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
		for (typename MeshVertexRange::Iterator vert = range.begin(); vert;
				++vert) {
			OperatorT::apply(*vert, mGrid, mDt);
		}
	}

protected:
	InterruptT* mInterrupt;

};
template<typename OperatorT,
		typename InterruptT = openvdb::util::NullInterrupter>
class ComputePertubationOperator {
public:
	SpringLevelSet& mGrid;
	ComputePertubationOperator(SpringLevelSet& grid, InterruptT* _interrupt,
			double t = 0.0, double dt = 0.0,
			TemporalIntegrationScheme scheme =
					TemporalIntegrationScheme::RK1) :
			mGrid(grid), mInterrupt(_interrupt), mIntegrationScheme(scheme), mTime(
					t), mTimeStep(dt) {
	}
	virtual ~ComputePertubationOperator() {
	}
	void process(bool threaded = true) {
		if (mInterrupt)
			mInterrupt->start("Processing springls");
		OperatorT::init(mGrid);
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
		for (typename SpringlRange::Iterator springl = range.begin(); springl; ++springl) {
			OperatorT::compute(*springl, mGrid, mTime);
		}
	}

protected:
	InterruptT* mInterrupt;
	TemporalIntegrationScheme mIntegrationScheme;
	double mTime;
	double mTimeStep;

};



#endif /* VDBTOOLS_INCLUDE_PERTURBATIONOPERATION_H_ */
