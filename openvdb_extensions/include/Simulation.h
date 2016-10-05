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

#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "SpringLevelSet.h"
#include <thread>
#include <mutex>
#include <chrono>
#include <AlloyParameterPane.h>
#include <AlloyWorker.h>
class Simulation;
class SimulationListener{
public:
	virtual void SimulationEvent(Simulation* simulation,int mSimulationIteration,double time)=0;
	virtual ~SimulationListener();
};
class Simulation: public aly::RecurrentTask {
protected:
	bool mPaused;
	double mComputeTimeSeconds;
	std::string mName;
	bool mIsInitialized;
	double mTimeStep;
	double mSimulationDuration;
	double mSimulationTime;
	std::string outputDirectory;
	uint64_t mSimulationIteration;
	std::thread mSimulationThread;
	virtual bool stepInternal()=0;
public:
	std::function<void(uint64_t iteration,bool lastIteration)> onUpdate;
	typedef std::chrono::high_resolution_clock Clock;
	virtual bool init()=0;
	virtual void cleanup()=0;
	virtual void setup(const aly::ParameterPanePtr& pane)=0;
	Simulation(const std::string& name);
	inline const std::string& getName()const {return mName;}
	inline void setName(const std::string& name){mName=name;}
	inline double getSimulationTime()const {return mSimulationTime;}
	inline double getSimulationDuration() const {return mSimulationDuration;}
	inline long getSimulationIteration()const {return mSimulationIteration;}
	inline double getComputeTimePerFrame() const {return mComputeTimeSeconds;}
	inline double getTimeStep() const {return mTimeStep;}
	virtual ~Simulation(){};
};

#endif /* SIMULATION_H_ */
