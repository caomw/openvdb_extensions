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
#ifndef INCLUDE_VIEWER_H_
#define INCLUDE_VIEWER_H_

#include <AlloyApplication.h>
#include <AlloyWidget.h>
#include <AlloyColorSelector.h>
#include <AlloyParameterPane.h>
#include <AlloyNumber.h>
#include <CommonShaders.h>
#include <GLFrameBuffer.h>
#include <memory>
#include <openvdb/openvdb.h>
class Viewer: public aly::Application {
protected:
	aly::Mesh mesh;
	aly::Mesh bboxMesh;
	aly::ImageRGBA renderImage;
	aly::Camera camera;
	std::string filePath;
	aly::CompositePtr renderRegion;
	int selectedIndex;
	int displayIndex;
	int cameraType;
	int shadingType;
	aly::TextLabelPtr objectNameLabel, vertexLabel, triangleLabel, quadLabel, colorLabel, texLabel;
	aly::Number lineWidth;
	aly::Number particleSize;
	aly::Color surfaceColor;
	aly::Color faceColor;
	aly::Color lineColor;
	aly::Color pointColor;
	std::string matcapImageFile;
	bool frameBuffersDirty;
	bool parametersDirty;
	aly::Number minDepth;
	aly::Number maxDepth;
	aly::int2 lastDepth;
	openvdb::FloatGrid::Ptr vdbGrid;
	aly::SelectionPtr displayIndexField;
	aly::SelectionPtr shadingStyleField;
	aly::ModifiableNumberPtr lineWidthField;
	aly::ModifiableNumberPtr particleSizeField;
	std::pair<aly::ModifiableNumberPtr, aly::ModifiableNumberPtr> rangerSliderField;
	aly::ColorSelectorPtr surfaceColorField;
	aly::ColorSelectorPtr pointColorField;
	aly::ColorSelectorPtr faceColorField;
	aly::ColorSelectorPtr lineColorField;

	std::unique_ptr<aly::GLFrameBuffer> colorFrameBuffer;
	std::unique_ptr<aly::GLFrameBuffer> depthFrameBuffer;
	std::unique_ptr<aly::GLFrameBuffer> compositeBuffer;
	std::unique_ptr<aly::GLFrameBuffer> wireframeFrameBuffer;

	std::unique_ptr<aly::GLFrameBuffer> lineFrameBuffer;
	std::unique_ptr<aly::CompositeShader> compositeShader;
	std::unique_ptr<aly::DepthAndNormalShader> depthAndNormalShader;
	std::unique_ptr<aly::ColorVertexShader> colorVertexShader;
	std::unique_ptr<aly::ParticleDepthShader> particleDepthShader;
	std::unique_ptr<aly::ParticleMatcapShader> particleMatcapShader;
	std::unique_ptr<aly::WireframeShader> wireframeShader;
	std::unique_ptr<aly::LineDistanceShader> lineDistanceShader;
	std::unique_ptr<aly::MatcapShader> matcapShader;
	std::unique_ptr<aly::ImageShader> imageShader;

	void initializeFrameBuffers(aly::AlloyContext* context);

	void load();
public:
	static const aly::box3f renderBBox;
	Viewer();
	void open(const std::string& filePath) {
		this->filePath = filePath;
	}
	bool init(aly::Composite& rootNode);
	void draw(aly::AlloyContext* context);
};

#endif /* INCLUDE_VIEWER_H_ */
