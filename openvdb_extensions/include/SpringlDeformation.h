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

#ifndef INCLUDE_SPRINGLDEFORMATION_H_
#define INCLUDE_SPRINGLDEFORMATION_H_
#include <openvdb/openvdb.h>
#include <openvdb/tools/LevelSetAdvect.h>
#include <openvdb/tools/LevelSetTracker.h>

#include <AlloyMath.h>
#include "SpringLevelSet.h"
#include "SpringlRange.h"
#include "MaxVelocityOperation.h"
template<typename ParticleAdvectionFunc,typename InterruptT = openvdb::util::NullInterrupter>
class SpringLevelSetParticleDeformation {
private:
	ParticleAdvectionFunc& mParticleAdvection;
	SpringLevelSet& mGrid;
	bool mResample;
	int mSignChanges;
	float mConvergenceThresold=0.01;
	int mTrackingIterations=128;
	double dt = 0.0;
	double time=0.0;
	std::mutex mSignChangeLock;
public:
	typedef openvdb::FloatGrid GridType;
	typedef openvdb::tools::LevelSetTracker<openvdb::FloatGrid, InterruptT> TrackerT;
	typedef typename TrackerT::RangeType RangeType;
	typedef typename TrackerT::LeafType LeafType;
	typedef typename TrackerT::BufferType BufferType;
	typedef typename TrackerT::ValueType ScalarType;
	void setConvergenceThreshold(float convg){
		mConvergenceThresold=convg;
	}
	void setTrackingIterations(int iters){
		mTrackingIterations=iters;
	}
	TemporalIntegrationScheme mTemporalScheme;
	MotionScheme mMotionScheme;
	InterruptT* mInterrupt;
	// disallow copy by assignment
	void operator=(const SpringLevelSetParticleDeformation& other) {
	}
	SpringLevelSetParticleDeformation(SpringLevelSet& grid,ParticleAdvectionFunc& advectFunc,
			MotionScheme scheme =
					MotionScheme::SEMI_IMPLICIT,
			InterruptT* interrupt = NULL) :mParticleAdvection(advectFunc),
			mSignChanges(0), mMotionScheme(scheme), mGrid(grid), mInterrupt(interrupt), mTemporalScheme(
					TemporalIntegrationScheme::RK4b), mResample(true) {
		mGrid.mConstellation.mParticleVelocity.resize(mGrid.mConstellation.mParticles.size(),openvdb::Vec3s(0.0));
		mGrid.mConstellation.mVertexVelocity.resize(mGrid.mConstellation.mVertexes.size(),openvdb::Vec3s(0.0));
		mGrid.mConstellation.mParticleLabel.resize(mGrid.mConstellation.mParticles.size(),0);
	}
	/// @return the temporal integration scheme
	TemporalIntegrationScheme getTemporalScheme() const {
		return mTemporalScheme;
	}
	/// @brief Set the spatial finite difference scheme
	void setTemporalScheme(TemporalIntegrationScheme scheme) {
		mTemporalScheme = scheme;
	}
	/// @brief Set enable resampling
	void setResampleEnabled(bool resample) {
		mResample = resample;
	}

	void advect(double startTime, double endTime) {
			const openvdb::math::Transform& trans = mGrid.mSignedLevelSet->transform();
			if (trans.mapType() == openvdb::math::UniformScaleMap::mapType()) {
				advect1<openvdb::math::UniformScaleMap>(startTime, endTime);
			} else if (trans.mapType()== openvdb::math::UniformScaleTranslateMap::mapType()) {
				advect1<openvdb::math::UniformScaleTranslateMap>(startTime,endTime);
			} else if (trans.mapType() == openvdb::math::UnitaryMap::mapType()) {
				advect1<openvdb::math::UnitaryMap>(startTime, endTime);
			} else if (trans.mapType() == openvdb::math::TranslationMap::mapType()) {
				advect1<openvdb::math::TranslationMap>(startTime, endTime);
			}
	}
	template<typename MapT> void track(double time) {
		const int RELAX_OUTER_ITERS = 1;
		const int RELAX_INNER_ITERS = 5;
		mGrid.updateUnSignedLevelSet();
		mGrid.updateNearestNeighbors();
		//Need this for original method
		//for (int iter = 0; iter < RELAX_OUTER_ITERS; iter++) {
			//mGrid.updateNearestNeighbors();
			//mGrid.relax(RELAX_INNER_ITERS);
		//}
		static int counter=0;
		if (mMotionScheme == MotionScheme::SEMI_IMPLICIT) {
			mGrid.updateUnSignedLevelSet(2.5 * openvdb::LEVEL_SET_HALF_WIDTH);
			mGrid.updateGradient();
			//WriteToRawFile(mGrid.mUnsignedLevelSet,MakeString()<<"/home/blake/tmp/unsigned"<<counter);
			TrackerT mTracker(*mGrid.mSignedLevelSet, mInterrupt);
			SpringLevelSetEvolve<MapT> evolve(*this, mTracker, time, 0.75, mTrackingIterations,mConvergenceThresold);
			evolve.process();
			//WriteToRawFile(mGrid.mSignedLevelSet,MakeString()<<"/home/blake/tmp/signed_after"<<counter);
			counter++;
		} else if (mMotionScheme == MotionScheme::EXPLICIT) {
			mGrid.mIsoSurface.updateVertexNormals(0);
			mGrid.mIsoSurface.dilate(0.5f);
			mGrid.updateSignedLevelSet();
			mGrid.updateUnSignedLevelSet(2.5 * openvdb::LEVEL_SET_HALF_WIDTH);
			mGrid.updateGradient();
			TrackerT mTracker(*mGrid.mSignedLevelSet, mInterrupt);
			SpringLevelSetEvolve<MapT> evolve(*this, mTracker, time, 0.75, mTrackingIterations,mConvergenceThresold);
			evolve.process();
		}
		if (mResample) {
			int cleaned = mGrid.clean();
			mGrid.updateUnSignedLevelSet();
			mGrid.updateIsoSurface();
			int added=mGrid.fill();
			mGrid.fillWithNearestNeighbors();
			//std::cout<<"Particle Filled "<<added<<" "<<100*added/(double)mGrid.mConstellation.getNumSpringls()<<"% "<<std::endl;
		} else {
			mGrid.updateIsoSurface();
		}
	}
	template<typename MapT> void evolve1(){
		TrackerT mTracker(*mGrid.mSignedLevelSet, mInterrupt);
		SpringLevelSetEvolve<MapT> ev(*this, mTracker, 0, 0.75, mTrackingIterations,mConvergenceThresold);
		ev.process();
	}
	template<typename MapT> void advect1(double mStartTime, double mEndTime,bool threaded=true) {
		openvdb::Vec3d vsz = mGrid.transformPtr()->voxelSize();
		double scale = std::max(std::max(vsz[0], vsz[1]), vsz[2]);
		const double EPS = 1E-30f;
		double voxelDistance = 0;
		const double MAX_TIME_STEP = SpringLevelSet::MAX_VEXT;
		mGrid.resetMetrics();
		for (time = mStartTime; time < mEndTime; time += dt) {
			MaxParticleVelocityOperator<InterruptT> op2(mGrid.mConstellation,mInterrupt);
			double err=std::sqrt(op2.process());
			double maxV = std::max(EPS, err);
			dt = aly::clamp(MAX_TIME_STEP * scale / std::max(1E-30, maxV), 0.0,mEndTime - time);
			if (dt < EPS) {
				break;
			}
			int N=mGrid.mConstellation.springls.size();
			openvdb::math::Transform::Ptr trans = mGrid.transformPtr();
			SpringlRange range(mGrid.mConstellation);
			if(threaded){
				tbb::parallel_for(range,*this);
			} else {
				(*this)(range);
			}
			//need this! commented out for debugging.
			if (mMotionScheme == MotionScheme::SEMI_IMPLICIT)track<MapT>(time);
		}

		if (mMotionScheme == MotionScheme::EXPLICIT)track<MapT>(time);
		mGrid.mConstellation.updateVertexNormals(0,0);
	}
	void operator()(const SpringlRange& range) const {
		if (openvdb::util::wasInterrupted(mInterrupt))
			tbb::task::self().cancel_group_execution();
		for (typename SpringlRange::Iterator springl = range.begin(); springl;++springl) {
			mParticleAdvection(*springl,time,dt);
		}
	}
	void evolve(){
		const openvdb::math::Transform& trans = mGrid.mSignedLevelSet->transform();
		if (trans.mapType() == openvdb::math::UniformScaleMap::mapType()) {
			evolve1<openvdb::math::UniformScaleMap>();
		} else if (trans.mapType()== openvdb::math::UniformScaleTranslateMap::mapType()) {
			evolve1<openvdb::math::UniformScaleTranslateMap>();
		} else if (trans.mapType() == openvdb::math::UnitaryMap::mapType()) {
			evolve1<openvdb::math::UnitaryMap>();
		} else if (trans.mapType() == openvdb::math::TranslationMap::mapType()) {
			evolve1<openvdb::math::TranslationMap>();
		}
	}
	template<typename MapT> class SpringLevelSetEvolve {
	public:
		SpringLevelSetParticleDeformation& mParent;
		typename TrackerT::LeafManagerType& mLeafs;
		TrackerT& mTracker;
		openvdb::tools::DiscreteField<openvdb::VectorGrid> mDiscreteField;
		const MapT* mMap;
		ScalarType mDt;
		double mTime;
		double mTolerance;
		int mIterations;
		SpringLevelSetEvolve(SpringLevelSetParticleDeformation& parent, TrackerT& tracker,
				double time, double dt, int iterations, double tolerance) :
				mMap(NULL), mParent(parent), mTracker(tracker), mIterations(
						iterations), mDiscreteField(*parent.mGrid.mGradient), mTime(
						time), mDt(dt), mTolerance(tolerance), mLeafs(
						tracker.leafs()) {
			mParent.mSignChanges = 0;
		}
		void process(bool threaded = true) {
			mMap = (mTracker.grid().transform().template constMap<MapT>().get());
			if (mParent.mInterrupt)
			mParent.mInterrupt->start("Processing voxels");
			mParent.mSignChanges=0;
			const int MIN_NUM_SIGN_CHANGES=32;
			int maxSignChanges=MIN_NUM_SIGN_CHANGES;
			int iter;
			for(iter=0;iter<mIterations;iter++) {
				mLeafs.rebuildAuxBuffers(1);
				mParent.mSignChanges=0;
				if (threaded) {
					tbb::parallel_for(mLeafs.getRange(mTracker.getGrainSize()), *this);
				} else {
					(*this)(mLeafs.getRange(mTracker.getGrainSize()));
				}
				mLeafs.swapLeafBuffer(1, mTracker.getGrainSize()==0);
				mLeafs.removeAuxBuffers();
				mTracker.track();
				maxSignChanges=std::max(mParent.mSignChanges,maxSignChanges);
				float ratio=(mParent.mSignChanges/(float)maxSignChanges);
				if(ratio<mTolerance){
					break;
				}
			}
			if (mParent.mInterrupt){
				mParent.mInterrupt->end();
			}
		}
		void operator()(const typename TrackerT::LeafManagerType::RangeType& range) const {
			using namespace openvdb;
			typedef math::BIAS_SCHEME<math::BiasedGradientScheme::FIRST_BIAS> Scheme;
			typedef typename Scheme::template ISStencil<FloatGrid>::StencilType Stencil;
			const math::Transform& trans = mTracker.grid().transform();
			typedef typename LeafType::ValueOnCIter VoxelIterT;
			const MapT& map = *mMap;
			Stencil stencil(mTracker.grid());
			int count=0;
			int signChanges=0;
			for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
				BufferType& result = mLeafs.getBuffer(n, 1);
				for (VoxelIterT iter = mLeafs.leaf(n).cbeginValueOn(); iter;++iter) {
					stencil.moveTo(iter);
					const Vec3s V = mDiscreteField(map.applyMap(iter.getCoord().asVec3d()), mTime);
					const Vec3s G = math::GradientBiased<MapT,openvdb::math::BiasedGradientScheme::FIRST_BIAS>::result(map, stencil, V);
					ScalarType delta=mDt * V.dot(G);
					ScalarType old=*iter;
					//Number of sign changes is a good indicator of the interface is moving.
					if(old*(old-delta)<0){
						signChanges++;
					}
					result.setValue(iter.pos(), old -  delta);
				}
			}
			mParent.mSignChangeLock.lock();
				mParent.mSignChanges+=signChanges;
			mParent.mSignChangeLock.unlock();
		}
	};
};

#endif /* INCLUDE_SPRINGLDEFORMATION_H_ */
