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
#include "Viewer.h"
#include "UtilitiesVDB.h"
#include <openvdb/openvdb.h>
using namespace aly;
using namespace openvdb;
const box3f Viewer::renderBBox = box3f(float3(-0.5f, -0.5f, -0.5f), float3(1.0f, 1.0f, 1.0f));
Viewer::Viewer() :
		Application(1200, 600, "OpenVDB Viewer"),selectedIndex(-1),displayIndex(0),cameraType(0),shadingType(0),frameBuffersDirty(true),parametersDirty(true){
}
void Viewer::initializeFrameBuffers(aly::AlloyContext* context) {
	float2 dims = renderRegion->getBounds().dimensions;
	int w = (int) dims.x;
	int h = (int) dims.y;
	compositeBuffer->initialize(w,h);
	colorFrameBuffer->initialize(w, h);
	depthFrameBuffer->initialize(w, h);
	wireframeFrameBuffer->initialize(w, h);
	lineFrameBuffer->initialize(w,h);
}
bool Viewer::init(Composite& rootNode) {
	openvdb::initialize();
	lastDepth=int2(-1,-1);
	displayIndex = 0;
	cameraType=0;
	shadingType=0;
	parametersDirty = true;
	frameBuffersDirty = true;
	compositeBuffer.reset(new GLFrameBuffer());
	depthFrameBuffer.reset(new GLFrameBuffer());
	wireframeFrameBuffer.reset(new GLFrameBuffer());
	lineFrameBuffer.reset(new GLFrameBuffer());
	colorFrameBuffer.reset(new GLFrameBuffer());
	depthAndNormalShader.reset(new DepthAndNormalShader());
	wireframeShader.reset(new WireframeShader());
	particleDepthShader.reset(new ParticleDepthShader());
	particleMatcapShader.reset(new ParticleMatcapShader());
	matcapShader.reset(new MatcapShader());
	imageShader.reset(new ImageShader());
	lineDistanceShader.reset(new LineDistanceShader());
	colorVertexShader.reset(new ColorVertexShader());
	compositeShader.reset(new CompositeShader());

	matcapImageFile = getFullPath("images/JG_Silver.png");

	BorderCompositePtr layout = BorderCompositePtr(new BorderComposite("UI Layout", CoordPX(0.0f, 0.0f), CoordPercent(1.0f, 1.0f), false));
	ParameterPanePtr controls = ParameterPanePtr(new ParameterPane("Controls", CoordPX(0.0f, 0.0f), CoordPercent(1.0f, 1.0f)));
	BorderCompositePtr controlLayout = BorderCompositePtr(new BorderComposite("Control Layout", CoordPX(0.0f, 0.0f), CoordPercent(1.0f, 1.0f), true));

	controls->onChange = [this](const std::string& label,const AnyInterface& value) {
		if(label=="Matcap") {
			AlloyApplicationContext()->addDeferredTask([this] {
						matcapShader->setTextureImage(matcapImageFile);
					});
		}
		if(label=="Tree Depth"){
			if(vdbGrid.get()!=nullptr&&(lastDepth.x!=minDepth.toInteger()||lastDepth.y!=maxDepth.toInteger())){
				lastDepth.x=minDepth.toInteger();
				lastDepth.y=maxDepth.toInteger();
				AlloyApplicationContext()->addDeferredTask([this] {
					ConvertLevelSetToBBoxTree(vdbGrid,bboxMesh,minDepth.toInteger(),maxDepth.toInteger());
					camera.setPose(MakeTransform(bboxMesh.getBoundingBox(), renderBBox));
				});
			}

		}
		if(label=="Data") {
			load();
		}
		parametersDirty=true;
	};
	float aspect = 6.0f;
	lineWidth = Float(2.0f);
	particleSize = Float(0.2f);

	lineColor = AlloyApplicationContext()->theme.DARK.toSemiTransparent(0.5f);
	faceColor = AlloyApplicationContext()->theme.LIGHT;
	surfaceColor = Color(255, 255, 255, 255);
	pointColor = Color(255, 255, 255, 255);

	controls->setAlwaysShowVerticalScrollBar(false);
	controls->setScrollEnabled(false);
	controls->backgroundColor = MakeColor(getContext()->theme.DARKER);
	controlLayout->backgroundColor = MakeColor(getContext()->theme.DARKER);
	controlLayout->borderWidth = UnitPX(1.0f);
	controlLayout->borderColor = MakeColor(getContext()->theme.LIGHT);
	renderRegion = CompositePtr(new Composite("View", CoordPX(0.0f, 0.0f), CoordPercent(1.0f, 1.0f)));
	layout->setWest(controlLayout, UnitPX(300.0f));
	controlLayout->setNorth(controls, UnitPX(500.0f));
	layout->setCenter(renderRegion);
	CompositePtr infoComposite = CompositePtr(new Composite("Info", CoordPX(0.0f, 0.0f), CoordPercent(1.0f, 1.0f)));
	infoComposite->setOrientation(Orientation::Vertical);
	infoComposite->backgroundColor = MakeColor(getContext()->theme.DARKER);
	controlLayout->setSouth(infoComposite, UnitPX(170.0f));
	rootNode.add(layout);
	matcapShader->setTextureImage(matcapImageFile);
	ImageRGBA tmpImg(4, 4);
	tmpImg.set(RGBA(255, 255, 255, 255));
	particleMatcapShader->setTextureImage(tmpImg);
	camera.setNearFarPlanes(-10.0f, 10.0f);
	camera.setZoom(0.75f);
	camera.setCameraType(CameraType::Orthographic);

	auto dataField=controls->addFileField("Data", filePath, aspect);
	dataField->addFileExtensionRule("OpenVDB","vdb");
	dataField->addFileExtensionRule("Meshlab PLY","ply");
	dataField->addFileExtensionRule("Wavefront","obj");
	dataField->setFileExtensionRule(1);
	auto matcapFileField=controls->addFileField("Matcap", matcapImageFile, aspect);
	matcapFileField->addFileExtensionRule("Portable Network Graphics","png");
	matcapFileField->addFileExtensionRule("JPEG","jpg");
	matcapFileField->addFileExtensionRule("OpenEXR", "exr");
	matcapFileField->setFileExtensionRule(1);
	displayIndexField = controls->addSelectionField("Display", displayIndex, std::vector<std::string> { "Solid",  "Solid & BBox Tree", "Solid & Wireframe", "Points & Wireframe",
			"Wireframe", "Points" }, aspect);
	shadingStyleField = controls->addSelectionField("Shading", shadingType, std::vector<std::string> { "Matcap", "Vertex Color"}, aspect);
	controls->addSelectionField("Camera", cameraType, std::vector<std::string> { "Orthographic", "Perspective" }, aspect);
	lineWidthField = controls->addNumberField("Line Width", lineWidth, Float(1.0f), Float(10.0f), 5.5f);
	minDepth=Integer(1);
	maxDepth=Integer(1);
	rangerSliderField = controls->addRangeField("Tree Depth",minDepth,maxDepth,Integer(0),Integer(3),6.5f);
	particleSizeField = controls->addNumberField("Particle Size", particleSize, Float(0.0f), Float(1.0f), 5.5f);
	surfaceColorField = controls->addColorField("Surface", surfaceColor);
	pointColorField = controls->addColorField("Point", pointColor);
	faceColorField = controls->addColorField("Face", faceColor);
	lineColorField = controls->addColorField("Line", lineColor);
	addListener(&camera);
	renderRegion->onPack = [this]() {
		camera.setDirty(true);
		frameBuffersDirty=true;
	};
	camera.setActiveRegion(renderRegion.get());
	wireframeShader->setFaceColor(Color(0.1f, 0.1f, 1.0f, 0.5f));
	wireframeShader->setEdgeColor(Color(1.0f, 0.8f, 0.1f, 1.0f));
	wireframeShader->setLineWidth(lineWidth.toFloat());
	lineDistanceShader->setSolid(false);
	if(filePath.size()>0){
		load();
	}
	return true;
}
void Viewer::load() {

	std::string ext=GetFileExtension(filePath);
	if(ext=="vdb"){
		openvdb::io::File file(filePath);
		file.open();
		openvdb::GridPtrVecPtr grids = file.getGrids();
		if (grids->empty()) {
			return;
		}
		GridBase::Ptr grid = grids->front();
		vdbGrid = boost::dynamic_pointer_cast<FloatGrid>(grid);
		ConvertLevelSetToMesh(vdbGrid, mesh);
		ConvertLevelSetToBBoxTree(vdbGrid,bboxMesh,minDepth.toInteger(),maxDepth.toInteger());
		camera.setPose(MakeTransform(bboxMesh.getBoundingBox(), renderBBox));
	} else if(ext=="ply"||ext=="obj"){
		ReadMeshFromFile(filePath,mesh);
		camera.setPose(MakeTransform(mesh.getBoundingBox(), renderBBox));
	}
	camera.setDirty(true);
}
void Viewer::draw(AlloyContext* context) {
	if (frameBuffersDirty) {
		initializeFrameBuffers(context);
		frameBuffersDirty = false;
	}
	if (parametersDirty) {
		if (cameraType == 0) {
			camera.setNearFarPlanes(-10.0f, 10.0f);
			camera.setCameraType(CameraType::Orthographic);
		} else {
			camera.setNearFarPlanes(0.001f, 10.0f);
			camera.setCameraType(CameraType::Perspective);
		}
	}
	std::list<std::pair<const aly::Mesh*, aly::float4x4>> drawList;
	drawList.push_back(std::pair<const aly::Mesh*, aly::float4x4>(&mesh, mesh.pose));
	float psize = 4*(0.0025f + particleSize.toFloat() * 0.05f);
	if (camera.isDirty() || parametersDirty) {
		wireframeShader->setLineWidth(lineWidth.toFloat());
		switch (displayIndex) {
		case 0:
			if (shadingType == 1) {
				colorVertexShader->draw(drawList, camera, *colorFrameBuffer);
			} else {
				depthAndNormalShader->draw(drawList, camera, *depthFrameBuffer);
			}
			wireframeShader->setSolid(false);
			break;
		case 1:

			if (shadingType == 1) {
				colorVertexShader->draw(drawList, camera, *colorFrameBuffer);
			} else {
				depthAndNormalShader->draw(drawList, camera, *depthFrameBuffer);
			}
			lineDistanceShader->setLineWidth(lineWidth.toFloat());
			lineDistanceShader->draw({std::pair<const aly::Mesh*, aly::float4x4>(&bboxMesh, mesh.pose)}, camera, *lineFrameBuffer);
			wireframeShader->setSolid(false);
			wireframeShader->setFaceColor(Color(0, 0, 0, 0));
			wireframeShader->setEdgeColor(lineColor);
			wireframeShader->draw({std::pair<const aly::Mesh*, aly::float4x4>(&bboxMesh, mesh.pose)}, camera, *wireframeFrameBuffer);
			compositeBuffer->begin();
			if (shadingType == 1) {
				imageShader->draw(colorFrameBuffer->getTexture(), compositeBuffer->getViewport(), 1.0f, false);
			} else {
				matcapShader->draw(depthFrameBuffer->getTexture(), camera,compositeBuffer->getViewport(),compositeBuffer->getViewport(), surfaceColor.toRGBAf());
			}
			compositeBuffer->end();
			break;
		case 2:
			if (shadingType == 1) {
				colorVertexShader->draw(drawList, camera, *colorFrameBuffer);
			} else {
				depthAndNormalShader->draw(drawList, camera, *depthFrameBuffer);
			}
			wireframeShader->setSolid(true);
			wireframeShader->setFaceColor(Color(0, 0, 0, 0));
			wireframeShader->setEdgeColor(lineColor);
			wireframeShader->draw(drawList, camera, *wireframeFrameBuffer);
			break;
		case 3:
			if (shadingType == 1) {
				colorFrameBuffer->begin();
				particleMatcapShader->draw(drawList, camera, colorFrameBuffer->getViewport(), colorFrameBuffer->getViewport(), psize);
				colorFrameBuffer->end();
			} else {
				particleDepthShader->draw(drawList, camera, *depthFrameBuffer, psize);
			}
			wireframeShader->setSolid(false);
			wireframeShader->setFaceColor(Color(0, 0, 0, 0));
			wireframeShader->setEdgeColor(lineColor);
			wireframeShader->draw(drawList, camera, *wireframeFrameBuffer);
			break;
		case 4:
			depthAndNormalShader->draw(drawList, camera, *depthFrameBuffer);
			colorVertexShader->draw(drawList, camera, *depthFrameBuffer);
			wireframeShader->setSolid(false);
			wireframeShader->setFaceColor(Color(0, 0, 0, 0));
			wireframeShader->setEdgeColor(lineColor);
			wireframeShader->draw(drawList, camera, *wireframeFrameBuffer);
			break;
		case 5:
			if (shadingType == 1) {
				colorFrameBuffer->begin();
				particleMatcapShader->draw(drawList, camera, colorFrameBuffer->getViewport(), colorFrameBuffer->getViewport(), psize);
				colorFrameBuffer->end();
			} else {
				particleDepthShader->draw(drawList, camera, *depthFrameBuffer, psize);
			}
			break;
		default:
			break;
		}
	}
	glEnable(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
	const RGBAf bgColor = context->theme.DARKEST.toRGBAf();
	glClearColor(bgColor.x, bgColor.y, bgColor.z, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	box2px rbbox=renderRegion->getBounds();
	switch (displayIndex) {
	case 0:
		if (shadingType == 1) {
			imageShader->draw(colorFrameBuffer->getTexture(), context->pixelRatio * rbbox, 1.0f, false);
		} else {
			matcapShader->draw(depthFrameBuffer->getTexture(), camera, context->pixelRatio * rbbox, context->getViewport(), surfaceColor.toRGBAf());
		}
		break;
	case 1:

		compositeShader->draw(compositeBuffer->getTexture(),depthFrameBuffer->getTexture(),wireframeFrameBuffer->getTexture(),lineFrameBuffer->getTexture(), context->pixelRatio * rbbox);
		break;
	case 2:
		if (shadingType == 1) {
			imageShader->draw(colorFrameBuffer->getTexture(), context->pixelRatio * rbbox, 1.0f, false);
		} else {
			matcapShader->draw(depthFrameBuffer->getTexture(), camera, context->pixelRatio * rbbox, context->getViewport(), faceColor.toRGBAf());
		}
		imageShader->draw(wireframeFrameBuffer->getTexture(), context->pixelRatio * rbbox, 1.0f, false);
		break;
	case 3:
		imageShader->draw(wireframeFrameBuffer->getTexture(), context->pixelRatio * rbbox, 1.0f, false);
		if (shadingType == 1) {
			imageShader->draw(colorFrameBuffer->getTexture(), context->pixelRatio * rbbox, 1.0f, false);
		} else {
			matcapShader->draw(depthFrameBuffer->getTexture(), camera, context->pixelRatio * rbbox, context->getViewport(), pointColor.toRGBAf());
		}
		break;
	case 4:
		imageShader->draw(wireframeFrameBuffer->getTexture(), context->pixelRatio * rbbox, 1.0f, false);
		break;
	case 5:
		if (shadingType == 1) {
			imageShader->draw(colorFrameBuffer->getTexture(), context->pixelRatio * rbbox, 1.0f, false);
		} else {
			matcapShader->draw(depthFrameBuffer->getTexture(), camera, context->pixelRatio * rbbox, context->getViewport(), pointColor.toRGBAf());
		}
		break;
	default:
		break;
	}
	camera.setDirty(false);
	parametersDirty = false;
}

