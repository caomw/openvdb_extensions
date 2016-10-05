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

#ifndef ENRIGHTSIMULATION_H_
#define ENRIGHTSIMULATION_H_
#include <openvdb/openvdb.h>
#include "Simulation.h"
#include "SpringlFieldDeformation.h"
#include "SpringlCache.h"
#include <memory>

class EnrightSimulation:public Simulation {
	typedef openvdb::tools::EnrightField<float> FieldT;
	typedef SpringLevelSetFieldDeformation<FieldT> AdvectT;
protected:
	FieldT mField;
	std::unique_ptr<AdvectT> mAdvect;
	aly::Number mGridSize;
	MotionScheme mMotionScheme;
	SpringLevelSet mSource;
	std::shared_ptr<SpringlCache> cache;
	virtual bool stepInternal() override;
public:
	EnrightSimulation(int gridSize=128,MotionScheme motionScheme=MotionScheme::SEMI_IMPLICIT,const	std::shared_ptr<SpringlCache>& cache=nullptr);
	virtual bool init()override;
	virtual void cleanup() override;
	virtual void setup(const aly::ParameterPanePtr& pane) override;
};

#endif /* ENRIGHTSIMULATION_H_ */
