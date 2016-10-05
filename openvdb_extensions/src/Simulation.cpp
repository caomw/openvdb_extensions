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

#include "Simulation.h"
#include <boost/filesystem.hpp>
#include <openvdb/openvdb.h>
#include <sstream>
#include <fstream>
#include <ostream>
using namespace openvdb;
using namespace openvdb::tools;
SimulationListener::~SimulationListener(){

}

Simulation::Simulation(const std::string& name):
RecurrentTask(
[this](uint64_t iteration){
	bool ret=stepInternal();
	if(onUpdate){
		onUpdate(iteration,!ret);
	}
	return ret;
},5),
		mPaused(false),
		mComputeTimeSeconds(0.0),
		mName(name),
		mIsInitialized(false),
		mTimeStep(0),
		mSimulationDuration(0),
		mSimulationTime(0),
		mSimulationIteration(0) {

}


