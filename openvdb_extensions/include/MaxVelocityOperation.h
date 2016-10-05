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

#ifndef VDBTOOLS_INCLUDE_MAXVELOCITYOPERATION_H_
#define VDBTOOLS_INCLUDE_MAXVELOCITYOPERATION_H_
#include <openvdb/openvdb.h>
template<typename OperatorT, typename FieldT, typename InterruptT = openvdb::util::NullInterrupter>
class MaxVelocityOperator {
public:
	double mMaxAbsV;
	SpringLevelSet& mGrid;
	MaxVelocityOperator(SpringLevelSet& grid, const FieldT& field, double t, InterruptT* _interrupt) :
			mGrid(grid), mField(field), mTime(t), mInterrupt(_interrupt), mMaxAbsV(std::numeric_limits<double>::min()) {

	}
	MaxVelocityOperator(MaxVelocityOperator& other, tbb::split) :
			mGrid(other.mGrid), mMaxAbsV(other.mMaxAbsV), mField(other.mField), mTime(other.mTime), mInterrupt(NULL) {
	}
	virtual ~MaxVelocityOperator() {
	}
	double process(bool threaded = true) {
		if (mInterrupt)
			mInterrupt->start("Processing springls");
		SpringlRange range(mGrid.mConstellation);
		if (threaded) {
			tbb::parallel_reduce(range, *this);
		} else {
			(*this)(range);
		}
		if (mInterrupt)
			mInterrupt->end();
		return mMaxAbsV;
	}
	void join(const MaxVelocityOperator& other) {
		mMaxAbsV = std::max(mMaxAbsV, other.mMaxAbsV);
	}

	/// @note Never call this public method directly - it is called by
	/// TBB threads only!
	void operator()(const SpringlRange& range) {
		if (openvdb::util::wasInterrupted(mInterrupt))
			tbb::task::self().cancel_group_execution();
		OperatorT OpT;
		for (typename SpringlRange::Iterator springl = range.begin(); springl; ++springl) {
			mMaxAbsV = std::max(mMaxAbsV, OpT.findTimeStep(*springl, mGrid, mField, mTime));
		}
	}

protected:
	double mTime;
	InterruptT* mInterrupt;
	const FieldT& mField;
};
template<typename InterruptT = openvdb::util::NullInterrupter>
class MaxParticleVelocityOperator {
public:
	double mMaxAbsV;
	Constellation& mConstellation;
	InterruptT* mInterrupt;
	MaxParticleVelocityOperator(Constellation& constellation, InterruptT* _interrupt) :
			mConstellation(constellation), mInterrupt(_interrupt), mMaxAbsV(std::numeric_limits<double>::min()) {
	}
	MaxParticleVelocityOperator(MaxParticleVelocityOperator& other, tbb::split) :
			mConstellation(other.mConstellation), mMaxAbsV(other.mMaxAbsV), mInterrupt(other.mInterrupt) {
	}
	virtual ~MaxParticleVelocityOperator() {
	}
	double process(bool threaded = true) {
		if (mInterrupt)
			mInterrupt->start("Processing springls");
		SpringlRange range(mConstellation);
		if (threaded) {
			tbb::parallel_reduce(range, *this);
		} else {
			(*this)(range);
		}
		if (mInterrupt)
			mInterrupt->end();
		return mMaxAbsV;
	}
	void join(const MaxParticleVelocityOperator& other) {
		mMaxAbsV = std::max(mMaxAbsV, other.mMaxAbsV);

	}
	void operator()(const SpringlRange& range) {
		if (openvdb::util::wasInterrupted(mInterrupt))
			tbb::task::self().cancel_group_execution();
		for (typename SpringlRange::Iterator springl = range.begin(); springl; ++springl) {
			mMaxAbsV = std::max((double) springl->particleVelocity().lengthSqr(), mMaxAbsV);
		}
	}
};
template<typename FieldT, typename InterruptT = openvdb::util::NullInterrupter>
class MaxLevelSetVelocityOperator {
public:
	double mMaxAbsV;
	openvdb::FloatGrid& mGrid;
	typedef typename openvdb::tree::LeafManager<openvdb::FloatGrid::TreeType> LeafManagerType;
	LeafManagerType mLeafs;
	InterruptT* mInterrupt;
	double mTime;
	const FieldT& mField;
	MaxLevelSetVelocityOperator(openvdb::FloatGrid& grid, const FieldT& field, double t, InterruptT* _interrupt) :
			mGrid(grid), mField(field), mTime(t), mInterrupt(_interrupt), mMaxAbsV(std::numeric_limits<double>::min()), mLeafs(LeafManagerType(grid.tree())) {
	}
	MaxLevelSetVelocityOperator(MaxLevelSetVelocityOperator& other, tbb::split) :
			mGrid(other.mGrid), mMaxAbsV(other.mMaxAbsV), mField(other.mField), mTime(other.mTime), mInterrupt(other.mInterrupt), mLeafs(
					LeafManagerType(other.mGrid.tree())) {
	}
	virtual ~MaxLevelSetVelocityOperator() {
	}
	double process(bool threaded = true) {
		if (mInterrupt)
			mInterrupt->start("Processing springls");
		if (threaded) {
			tbb::parallel_reduce(mLeafs.getRange(1), *this);
		} else {
			(*this)(mLeafs.getRange(1));
		}
		if (mInterrupt)
			mInterrupt->end();
		return mMaxAbsV;
	}
	void join(const MaxLevelSetVelocityOperator& other) {
		mMaxAbsV = std::max(mMaxAbsV, other.mMaxAbsV);

	}
	void operator()(const typename LeafManagerType::RangeType& range) {
		if (openvdb::util::wasInterrupted(mInterrupt))
			tbb::task::self().cancel_group_execution();
		typedef typename LeafManagerType::LeafType::ValueOnCIter VoxelIterT;
		for (size_t n = range.begin(); n != range.end(); ++n) {
			for (VoxelIterT iter = mLeafs.leaf(n).cbeginValueOn(); iter; ++iter) {
				openvdb::math::Transform::Ptr trans = mGrid.transformPtr();
				openvdb::Vec3d pt = trans->indexToWorld(iter.getCoord());
				const openvdb::Vec3s V = mField(pt, mTime);
				mMaxAbsV = std::max(mMaxAbsV, (double) V.lengthSqr());
			}
		}
	}
};

#endif /* VDBTOOLS_INCLUDE_MAXVELOCITYOPERATION_H_ */
