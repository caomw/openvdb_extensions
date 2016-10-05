/*
 * Copyright(C) 2014, Blake C. Lucas, Ph.D. (img.science@gmail.com)
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

#include "EnrightSimulation.h"
#include <openvdb/openvdb.h>
#include <openvdb/tools/LevelSetUtil.h>
#include <openvdb/tools/LevelSetSphere.h>
EnrightSimulation::EnrightSimulation(int gridSize,MotionScheme scheme,const	std::shared_ptr<SpringlCache>& cache):Simulation("Enright"),cache(cache),mGridSize(aly::Integer(gridSize)),mMotionScheme(scheme) {

}
bool EnrightSimulation::init(){
	const float radius = 0.15f;
	const openvdb::Vec3f center(0.35f, 0.35f, 0.35f);
	float voxelSize = 1.0f / (float) (mGridSize.toInteger() - 1);
	openvdb::FloatGrid::Ptr mSignedLevelSet =openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(radius,center, voxelSize);
	mSource.create(*mSignedLevelSet);
	//Important! re-normalize distance to be in voxel units.
    openvdb::BBoxd bbox=mSource.mIsoSurface.updateBoundingBox();
    mAdvect=std::unique_ptr<AdvectT>(new AdvectT(mSource,mField,mMotionScheme));
	mAdvect->setTemporalScheme(TemporalIntegrationScheme::RK4b);
	mSimulationDuration=3.0f;
	mSimulationIteration=0;
	mTimeStep=0.5*voxelSize;
	mSource.setCurrentIteration(mSimulationIteration);
	std::string outFile=aly::MakeString()<<outputDirectory<<ALY_PATH_SEPARATOR<<"simulation_"<<std::setw(5)<<std::setfill('0')<<(mSimulationIteration+1)<<".json";
	if(aly::FileExists(outFile)){
		mSource.read(outFile);
	} else {
		mSource.write(outFile);
	}
	if(cache.get()!=nullptr){
		cache->set(mSimulationIteration,mSource);
	}
	return true;
}
void EnrightSimulation::cleanup(){
	mAdvect.reset();
}
void EnrightSimulation::setup(const aly::ParameterPanePtr& pane){
	mMotionScheme=MotionScheme::SEMI_IMPLICIT;
	outputDirectory=aly::GetDesktopDirectory()+ALY_PATH_SEPARATOR+"springls";
	aly::MakeDirectory(outputDirectory);
	pane->addDirectoryField("Output",outputDirectory,6.0f);
	pane->addSelectionField("Method",(int&)mMotionScheme,std::vector<std::string>{"Implicit","Semi-Implicit","Explicit"},6.0f);
	pane->addNumberField("Grid Size",mGridSize,4.0f);

}
bool EnrightSimulation::stepInternal(){
	Clock::time_point t0 = Clock::now();
	std::string outFile=aly::MakeString()<<outputDirectory<<ALY_PATH_SEPARATOR<<"simulation_"<<std::setw(5)<<std::setfill('0')<<(mSimulationIteration+1)<<".json";
	if(aly::FileExists(outFile)){
		mSource.read(outFile);
		mSimulationIteration++;
	} else {
		mAdvect->advect(mSimulationTime,mSimulationTime+mTimeStep);
		Clock::time_point t1 = Clock::now();
		mComputeTimeSeconds= 1E-6*std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
		mSimulationIteration++;
		mSource.setCurrentIteration(mSimulationIteration);
		mSource.write(outFile);
	}
	if(cache.get()!=nullptr){
		cache->set(mSimulationIteration,mSource);
	}
	mSimulationTime=mTimeStep*mSimulationIteration;
	if(mSimulationTime<=mSimulationDuration&&!isCanceled()){
		return true;
	} else {
		mSimulationIteration--;
		mSimulationTime=mSimulationDuration;
		return false;
	}
}

