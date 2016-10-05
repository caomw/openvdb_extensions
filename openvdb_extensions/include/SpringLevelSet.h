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
#ifndef INCLUDE_SPRINGLEVELSET_H_
#define INCLUDE_SPRINGLEVELSET_H_
#include "Springl.h"
#include "ParticleVolume.h"
#include "Constellation.h"
#include "MACGrid.h"
#include "RegularGrid.h"
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/openvdb.h>
#include <AlloyMesh.h>
enum class TemporalIntegrationScheme {
	UNKNOWN_TIS = openvdb::math::TemporalIntegrationScheme::UNKNOWN_TIS,
	RK1 = openvdb::math::TemporalIntegrationScheme::TVD_RK1,
	RK2 = openvdb::math::TemporalIntegrationScheme::TVD_RK2,
	RK3 = openvdb::math::TemporalIntegrationScheme::TVD_RK3,
	RK4a,
	RK4b
};
enum class MotionScheme {
	UNDEFINED=-1,
	IMPLICIT=0,
	SEMI_IMPLICIT=1,
	EXPLICIT=2
};
typedef std::vector<std::list<SpringlNeighbor>> NearestNeighborMap;
typedef openvdb::FloatGrid::Ptr SLevelSetPtr;
typedef openvdb::VectorGrid::Ptr SGradientPtr;
typedef openvdb::Int32Grid::Ptr SIndexPtr;
class SpringLevelSet {
protected:
	int mFillCount;
	int mCleanCount;
	int64_t mIteration;
	openvdb::tools::VolumeToMesh mVolToMesh;
	openvdb::math::Transform::Ptr mTransform;
	std::list<int> fillList;
	std::string constellationFile;
	std::string isoSurfaceFile;
	std::string signedLevelSetFile;
	std::string unsignedLevelSetFile;

public:
	static const float NEAREST_NEIGHBOR_RANGE; //voxel units
	static const int MAX_NEAREST_NEIGHBORS;
	static const float PARTICLE_RADIUS;
	static const float MAX_VEXT;
	static const float FILL_DISTANCE;
	static const float CLEAN_DISTANCE;
	static const float SHARPNESS;
	static const float SPRING_CONSTANT;
	static const float RELAX_TIMESTEP;
	static const float MIN_ASPECT_RATIO;
	static const float MAX_AREA;
	static const float MIN_AREA;
	Constellation mIsoSurface;
	ParticleVolume mParticleVolume;
	Constellation mConstellation;
	NearestNeighborMap mNearestNeighbors;
	SLevelSetPtr mSignedLevelSet;
	SLevelSetPtr mUnsignedLevelSet;
	SGradientPtr mGradient;
	SIndexPtr mSpringlIndexGrid;
	inline openvdb::math::Transform& transform() {
		return *mTransform;
	}
	inline openvdb::math::Transform::Ptr transformPtr() {
		return mTransform;
	}
	void setCurrentIteration(int64_t iter){
		mIteration=iter;
	}
	int64_t getCurrentIteration() const {
		return mIteration;
	}
	std::string getConstellationFile() const {
		return constellationFile;
	}
	std::string getIsoSurfaceFile() const {
		return isoSurfaceFile;
	}
	Springl& getSpringl(const openvdb::Index32 id);
	openvdb::Vec3s& getParticle(const openvdb::Index32 id);
	openvdb::Vec3s& getParticleNormal(const openvdb::Index32 id);
	openvdb::Vec3s& getSpringlVertex(const openvdb::Index32 id, const int i);
	openvdb::Vec3s& getSpringlVertex(const openvdb::Index32 gid);
	std::list<SpringlNeighbor>& getNearestNeighbors(openvdb::Index32 id, int8_t e);
	inline int getLastFillCount() const {
		return mFillCount;
	}
	inline int getLastCleanCount() const {
		return mCleanCount;
	}
	void resetMetrics() {
		mCleanCount = 0;
		mFillCount = 0;
	}
	int clean();
	int fill();
	void fillWithNearestNeighbors();
	void fillWithVelocityField(MACGrid<float>& grid, float radius);
	void evolve();
	void updateLines();
	void updateGradient();
	void updateIsoSurface();
	void saveResources(const std::string& directory);
	void loadResources(const std::string& directory);
	template<class Archive> void save(Archive& ar) const {
		ar(		cereal::make_nvp("iteration",mIteration),
				cereal::make_nvp("constellation",constellationFile),
				cereal::make_nvp("iso-surface",isoSurfaceFile),
				cereal::make_nvp("signed_levelset",signedLevelSetFile),
				cereal::make_nvp("unsigned_levelset",unsignedLevelSetFile));
	}
	template<class Archive> void load(Archive& ar) {
		ar(		cereal::make_nvp("iteration",mIteration),
				cereal::make_nvp("constellation",constellationFile),
				cereal::make_nvp("iso-surface",isoSurfaceFile),
				cereal::make_nvp("signed_levelset",signedLevelSetFile),
				cereal::make_nvp("unsigned_levelset",unsignedLevelSetFile));
	}
	void write(const std::string& file);
	void read(const std::string& file);
	void updateUnSignedLevelSet(double distance = openvdb::LEVEL_SET_HALF_WIDTH);
	void updateSignedLevelSet();
	void computeStatistics(const Constellation& mesh, openvdb::FloatGrid& levelSet);
	void computeStatistics(const Constellation& mesh);
	void relax(int iters = 10);
	double distanceToConstellation(const openvdb::Vec3s& pt);
	void updateNearestNeighbors(bool threaded = true);
	void create(const Constellation& mesh, openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform());
	void create(openvdb::FloatGrid& grid);
	void create(RegularGrid<float>& grid);
	SpringLevelSet() :
		mFillCount(0),mCleanCount(0), mIteration(0), mVolToMesh(0.0), mTransform(openvdb::math::Transform::createLinearTransform(1.0)) {
	}

	~SpringLevelSet() {
	}
};

#endif /* INCLUDE_SPRINGLEVELSET_H_ */
