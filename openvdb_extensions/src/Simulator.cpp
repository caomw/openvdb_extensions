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
#include "Simulator.h"
#include "UtilitiesVDB.h"
#include "EnrightSimulation.h"
#include <openvdb/openvdb.h>
using namespace aly;
using namespace openvdb;
const box3f Simulator::renderBBox = box3f(float3(-0.5f, -0.5f, -0.5f), float3(1.0f, 1.0f, 1.0f));
Simulator::Simulator() :
		Application(1200, 600, "OpenVDB Simulator"), selectedIndex(-1), displayIndex(0), cameraType(0), shadingType(0), frameBuffersDirty(true), parametersDirty(
				true) {
}
void Simulator::initializeFrameBuffers(aly::AlloyContext* context) {
	float2 dims = renderRegion->getBounds().dimensions;
	int w = (int) dims.x;
	int h = (int) dims.y;
	compositeBuffer->initialize(w, h);
	colorFrameBuffer->initialize(w, h);
	depthFrameBuffer->initialize(w, h);
	elemDepthFrameBuffer->initialize(w,h);
	wireframeFrameBuffer->initialize(w, h);
	lineFrameBuffer->initialize(w, h);
}
bool Simulator::init(Composite& rootNode) {
	openvdb::initialize();
	cache=std::shared_ptr<SpringlCache>(new SpringlCache());
	simulation = std::shared_ptr<Simulation>(new EnrightSimulation(256,MotionScheme::SEMI_IMPLICIT,cache));
	simulation->onUpdate=[this](uint64_t iteration,bool lastIteration){
		if(lastIteration){
			stopButton->setVisible(false);
			playButton->setVisible(true);
		}
		AlloyApplicationContext()->addDeferredTask([this](){
			timelineSlider->setUpperValue((int)simulation->getSimulationIteration());
			timelineSlider->setTimeValue((int)simulation->getSimulationIteration());

		});
	};
	lastDepth = int2(-1, -1);
	displayIndex = 0;
	cameraType = 0;
	shadingType = 0;
	parametersDirty = true;
	frameBuffersDirty = true;
	compositeBuffer.reset(new GLFrameBuffer());
	depthFrameBuffer.reset(new GLFrameBuffer());
	elemDepthFrameBuffer.reset(new GLFrameBuffer());
	wireframeFrameBuffer.reset(new GLFrameBuffer());
	lineFrameBuffer.reset(new GLFrameBuffer());
	colorFrameBuffer.reset(new GLFrameBuffer());
	depthAndNormalShader.reset(new DepthAndNormalShader());
	wireframeShader.reset(new WireframeShader());
	particleDepthShader.reset(new ParticleDepthShader());
	particleMatcapShader.reset(new ParticleMatcapShader());
	matcapShader.reset(new MatcapShader());
	springlShader.reset(new SpringlShader());
	imageShader.reset(new ImageShader());
	lineDistanceShader.reset(new LineDistanceShader());
	colorVertexShader.reset(new ColorVertexShader());
	compositeShader.reset(new CompositeShader());

	matcapImageFile = getFullPath("images/JG_Silver.png");
	matcapElemImageFile = getFullPath("images/JG_Red.png");

	BorderCompositePtr layout = BorderCompositePtr(new BorderComposite("UI Layout", CoordPX(0.0f, 0.0f), CoordPercent(1.0f, 1.0f), false));
	ParameterPanePtr controls = ParameterPanePtr(new ParameterPane("Controls", CoordPX(0.0f, 0.0f), CoordPercent(1.0f, 1.0f)));
	BorderCompositePtr controlLayout = BorderCompositePtr(new BorderComposite("Control Layout", CoordPX(0.0f, 0.0f), CoordPercent(1.0f, 1.0f), true));

	controls->onChange = [this](const std::string& label,const AnyInterface& value) {
		if(label=="Foreground") {
			AlloyApplicationContext()->addDeferredTask([this] {
						matcapShader->setTextureImage(matcapImageFile);
						springlShader->setBackgroundTextureImage(matcapImageFile);

					});
		}
		if(label=="Background") {
					AlloyApplicationContext()->addDeferredTask([this] {
						springlShader->setForegroundTextureImage(matcapElemImageFile);
					});
		}
		if(label=="Tree Depth") {
			if(vdbGrid.get()!=nullptr&&(lastDepth.x!=minDepth.toInteger()||lastDepth.y!=maxDepth.toInteger())) {
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
	controls->borderColor = MakeColor(getContext()->theme.DARK);
	controls->borderWidth=UnitPX(1.0f);

	controlLayout->backgroundColor = MakeColor(getContext()->theme.DARKER);
	controlLayout->borderWidth = UnitPX(0.0f);
	renderRegion = CompositePtr(new Composite("View", CoordPX(0.0f, 0.0f), CoordPercent(1.0f, 1.0f)));
	layout->setWest(controlLayout, UnitPX(300.0f));
	controlLayout->setCenter(controls);
	layout->setCenter(renderRegion);
	CompositePtr infoComposite = CompositePtr(new Composite("Info", CoordPX(0.0f, 0.0f), CoordPercent(1.0f, 1.0f)));
	infoComposite->backgroundColor = MakeColor(getContext()->theme.DARKER);
	infoComposite->borderColor = MakeColor(getContext()->theme.DARK);
	infoComposite->borderWidth=UnitPX(0.0f);
	playButton=IconButtonPtr(new IconButton(0xf144,CoordPerPX(0.5f,0.5f,-35.0f,-35.0f),CoordPX(70.0f,70.0f)));
	stopButton=IconButtonPtr(new IconButton(0xf28d,CoordPerPX(0.5f,0.5f,-35.0f,-35.0f),CoordPX(70.0f,70.0f)));
	playButton->borderWidth=UnitPX(0.0f);
	stopButton->borderWidth=UnitPX(0.0f);
	playButton->backgroundColor=MakeColor(getContext()->theme.DARKER);
	stopButton->backgroundColor=MakeColor(getContext()->theme.DARKER);
	playButton->foregroundColor=MakeColor(0,0,0,0);
	stopButton->foregroundColor=MakeColor(0,0,0,0);
	playButton->iconColor=MakeColor(getContext()->theme.LIGHTER);
	stopButton->iconColor=MakeColor(getContext()->theme.LIGHTER);
	playButton->borderColor=MakeColor(getContext()->theme.LIGHTEST);
	stopButton->borderColor=MakeColor(getContext()->theme.LIGHTEST);
	playButton->onMouseDown=[this](AlloyContext* context,const InputEvent& e){
		if(e.button==GLFW_MOUSE_BUTTON_LEFT){
			stopButton->setVisible(true);
			playButton->setVisible(false);
			simulation->cancel();
			cache->clear();
			simulation->init();
			int maxIteration=(int)std::ceil(simulation->getSimulationDuration() / simulation->getTimeStep());
			timelineSlider->setTimeValue(0);
			timelineSlider->setMaxValue(maxIteration);
			timelineSlider->setVisible(true);
			simulation->execute();
			return true;
		}
		return false;
	};
	stopButton->onMouseDown=[this](AlloyContext* context,const InputEvent& e){
		if(e.button==GLFW_MOUSE_BUTTON_LEFT){
			stopButton->setVisible(false);
			playButton->setVisible(true);
			simulation->cancel();
			return true;
		}
		return false;
	};
	stopButton->setVisible(false);
	infoComposite->add(playButton);
	infoComposite->add(stopButton);
	controlLayout->setSouth(infoComposite, UnitPX(80.0f));
	rootNode.add(layout);
	matcapShader->setTextureImage(matcapImageFile);
	springlShader->setBackgroundTextureImage(matcapImageFile);
	springlShader->setForegroundTextureImage(matcapElemImageFile);

	ImageRGBA tmpImg(4, 4);
	tmpImg.set(RGBA(255, 255, 255, 255));
	particleMatcapShader->setTextureImage(tmpImg);
	camera.setNearFarPlanes(-10.0f, 10.0f);
	camera.setZoom(0.75f);
	camera.setCameraType(CameraType::Orthographic);
	controls->addGroup("Visualization",true);
	auto matcapFileField = controls->addFileField("Foreground", matcapImageFile, aspect);
	auto matcapElemFileField = controls->addFileField("Background", matcapElemImageFile, aspect);

	matcapFileField->addFileExtensionRule("Portable Network Graphics", "png");
	matcapFileField->addFileExtensionRule("JPEG", "jpg");
	matcapFileField->addFileExtensionRule("OpenEXR", "exr");
	matcapFileField->setFileExtensionRule(1);

	matcapElemFileField->addFileExtensionRule("Portable Network Graphics", "png");
	matcapElemFileField->addFileExtensionRule("JPEG", "jpg");
	matcapElemFileField->addFileExtensionRule("OpenEXR", "exr");
	matcapElemFileField->setFileExtensionRule(1);
	displayIndexField = controls->addSelectionField("Display", displayIndex, std::vector<std::string> { "Solid", "Solid & BBox Tree", "Solid & Wireframe",
			"Points & Wireframe", "Wireframe", "Points" }, aspect);
	shadingStyleField = controls->addSelectionField("Shading", shadingType, std::vector<std::string> { "Matcap", "Vertex Color" }, aspect);
	controls->addSelectionField("Camera", cameraType, std::vector<std::string> { "Orthographic", "Perspective" }, aspect);
	lineWidthField = controls->addNumberField("Line Width", lineWidth, Float(1.0f), Float(10.0f), 5.5f);
	minDepth = Integer(1);
	maxDepth = Integer(1);
	rangerSliderField = controls->addRangeField("Tree Depth", minDepth, maxDepth, Integer(0), Integer(3), 6.5f);
	particleSizeField = controls->addNumberField("Particle Size", particleSize, Float(0.0f), Float(1.0f), 5.5f);
	surfaceColorField = controls->addColorField("Surface", surfaceColor);
	pointColorField = controls->addColorField("Point", pointColor);
	faceColorField = controls->addColorField("Face", faceColor);
	lineColorField = controls->addColorField("Line", lineColor);

	controls->addGroup("Simulation",true);
	simulation->setup(controls);
	timelineSlider = TimelineSliderPtr(
			new TimelineSlider("Timeline", CoordPerPX(0.0f, 1.0f, 0.0f, -80.0f), CoordPerPX(1.0f, 0.0f, 0.0f, 80.0f), Integer(0),
					Integer(0), Integer(0)));
	timelineSlider->backgroundColor = MakeColor(AlloyApplicationContext()->theme.DARKER);
	timelineSlider->borderColor=MakeColor(AlloyApplicationContext()->theme.DARK);
	timelineSlider->borderWidth=UnitPX(0.0f);
	timelineSlider->onChangeEvent = [this](const Number& timeValue, const Number& lowerValue, const Number& upperValue) {

	};
	timelineSlider->setMajorTick(100);
	timelineSlider->setMinorTick(10);
	timelineSlider->setLowerValue(0);
	timelineSlider->setUpperValue(0);
	timelineSlider->setVisible(false);
	timelineSlider->setModifiable(false);
	renderRegion->add(timelineSlider);
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
	if (filePath.size() > 0) {
		load();
	}
	return true;
}
void Simulator::load() {

	std::string ext = GetFileExtension(filePath);
	if (ext == "vdb") {
		openvdb::io::File file(filePath);
		file.open();
		openvdb::GridPtrVecPtr grids = file.getGrids();
		if (grids->empty()) {
			return;
		}
		GridBase::Ptr grid = grids->front();
		vdbGrid = boost::dynamic_pointer_cast<FloatGrid>(grid);
		ConvertLevelSetToMesh(vdbGrid, mesh);
		ConvertLevelSetToBBoxTree(vdbGrid, bboxMesh, minDepth.toInteger(), maxDepth.toInteger());
		camera.setPose(MakeTransform(bboxMesh.getBoundingBox(), renderBBox));
	} else if (ext == "ply" || ext == "obj") {
		ReadMeshFromFile(filePath, mesh);
		camera.setPose(MakeTransform(mesh.getBoundingBox(), renderBBox));
	}
	camera.setDirty(true);
}
void Simulator::draw(AlloyContext* context) {
	int fr=timelineSlider->getTimeValue().toInteger();
	if(fr!=lastValidTime){
		std::shared_ptr<CacheElement> elem=cache->get(fr);
		if(elem.get()!=nullptr){
			elem->getConstellation()->copyTo(constellation);
			elem->getIsoSurface()->copyTo(mesh);
			float4x4 T=MakeTransform(mesh.getBoundingBox(),renderBBox);
			camera.setPose(T);
			lastValidTime=fr;
		}
	}
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
	std::list<std::pair<const aly::Mesh*, aly::float4x4>> drawList2;
	drawList.push_back(std::pair<const aly::Mesh*, aly::float4x4>(&mesh, mesh.pose));
	drawList2.push_back(std::pair<const aly::Mesh*, aly::float4x4>(&constellation, constellation.pose));
	float psize = 4 * (0.0025f + particleSize.toFloat() * 0.05f);
	if (camera.isDirty() || parametersDirty) {
		wireframeShader->setLineWidth(lineWidth.toFloat());
		depthAndNormalShader->draw(drawList2, camera, *elemDepthFrameBuffer);
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
			lineDistanceShader->draw( { std::pair<const aly::Mesh*, aly::float4x4>(&bboxMesh, mesh.pose) }, camera, *lineFrameBuffer);
			wireframeShader->setSolid(false);
			wireframeShader->setFaceColor(Color(0, 0, 0, 0));
			wireframeShader->setEdgeColor(lineColor);
			wireframeShader->draw( { std::pair<const aly::Mesh*, aly::float4x4>(&bboxMesh, mesh.pose) }, camera, *wireframeFrameBuffer);
			compositeBuffer->begin();
			if (shadingType == 1) {
				imageShader->draw(colorFrameBuffer->getTexture(), compositeBuffer->getViewport(), 1.0f, false);
			} else {
				if(constellation.vertexLocations.size()>0){
					springlShader->draw( elemDepthFrameBuffer->getTexture(),depthFrameBuffer->getTexture(),camera,
							compositeBuffer->getViewport(), compositeBuffer->getViewport(),
							surfaceColor.toRGBAf());
				} else {
				matcapShader->draw(depthFrameBuffer->getTexture(), camera, compositeBuffer->getViewport(), compositeBuffer->getViewport(),
						surfaceColor.toRGBAf());
				}
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
			if(constellation.vertexLocations.size()>0){
				wireframeShader->draw(drawList2, camera, *wireframeFrameBuffer);
			} else {
				wireframeShader->draw(drawList, camera, *wireframeFrameBuffer);
			}
			break;
		case 3:
			if (shadingType == 1) {
				colorFrameBuffer->begin();
				if(constellation.vertexLocations.size()>0){
					particleMatcapShader->draw(drawList2, camera, colorFrameBuffer->getViewport(), colorFrameBuffer->getViewport(), psize);
				} else {
					particleMatcapShader->draw(drawList, camera, colorFrameBuffer->getViewport(), colorFrameBuffer->getViewport(), psize);
				}
				colorFrameBuffer->end();
			} else {
				if(constellation.vertexLocations.size()>0){
					particleDepthShader->draw(drawList2, camera, *depthFrameBuffer, psize);
				} else {
					particleDepthShader->draw(drawList, camera, *depthFrameBuffer, psize);

				}
			}
			wireframeShader->setSolid(false);
			wireframeShader->setFaceColor(Color(0, 0, 0, 0));
			wireframeShader->setEdgeColor(lineColor);
			if(constellation.vertexLocations.size()>0){
				wireframeShader->draw(drawList2, camera, *wireframeFrameBuffer);
			} else {
				wireframeShader->draw(drawList, camera, *wireframeFrameBuffer);
			}
			break;
		case 4:
			depthAndNormalShader->draw(drawList, camera, *depthFrameBuffer);
			colorVertexShader->draw(drawList, camera, *depthFrameBuffer);
			wireframeShader->setSolid(false);
			wireframeShader->setFaceColor(Color(0, 0, 0, 0));
			wireframeShader->setEdgeColor(lineColor);
			if(constellation.vertexLocations.size()>0){
				wireframeShader->draw(drawList2, camera, *wireframeFrameBuffer);
			} else {
				wireframeShader->draw(drawList, camera, *wireframeFrameBuffer);
			}
			break;
		case 5:
			if (shadingType == 1) {
				colorFrameBuffer->begin();
				if(constellation.vertexLocations.size()>0){
					particleMatcapShader->draw(drawList2, camera, colorFrameBuffer->getViewport(), colorFrameBuffer->getViewport(), psize);
				} else {
					particleMatcapShader->draw(drawList, camera, colorFrameBuffer->getViewport(), colorFrameBuffer->getViewport(), psize);

				}
				colorFrameBuffer->end();
			} else {
				if(constellation.vertexLocations.size()>0){
					particleDepthShader->draw(drawList2, camera, *depthFrameBuffer, psize);
				} else {
					particleDepthShader->draw(drawList, camera, *depthFrameBuffer, psize);
				}
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
	box2px rbbox = renderRegion->getBounds();
	switch (displayIndex) {
	case 0:
		if (shadingType == 1) {
			imageShader->draw(colorFrameBuffer->getTexture(), context->pixelRatio * rbbox, 1.0f, false);
		} else {
			if(constellation.vertexLocations.size()>0){
				springlShader->draw( elemDepthFrameBuffer->getTexture(),depthFrameBuffer->getTexture(),camera,
						 context->pixelRatio * rbbox, context->getViewport(),
						surfaceColor.toRGBAf());
			} else {
				matcapShader->draw(depthFrameBuffer->getTexture(), camera, context->pixelRatio * rbbox, context->getViewport(), surfaceColor.toRGBAf());
			}
		}
		break;
	case 1:

		compositeShader->draw(compositeBuffer->getTexture(), depthFrameBuffer->getTexture(), wireframeFrameBuffer->getTexture(), lineFrameBuffer->getTexture(),
				context->pixelRatio * rbbox);
		break;
	case 2:
		if (shadingType == 1) {
			imageShader->draw(colorFrameBuffer->getTexture(), context->pixelRatio * rbbox, 1.0f, false);
		} else {
			if(constellation.vertexLocations.size()>0){
				springlShader->draw( elemDepthFrameBuffer->getTexture(),depthFrameBuffer->getTexture(),camera,
						 context->pixelRatio * rbbox, context->getViewport(),
						surfaceColor.toRGBAf());
			} else {
				matcapShader->draw(depthFrameBuffer->getTexture(), camera, context->pixelRatio * rbbox, context->getViewport(), surfaceColor.toRGBAf());
			}
		}
		imageShader->draw(wireframeFrameBuffer->getTexture(), context->pixelRatio * rbbox, 1.0f, false);
		break;
	case 3:
		imageShader->draw(wireframeFrameBuffer->getTexture(), context->pixelRatio * rbbox, 1.0f, false);
		if (shadingType == 1) {
			imageShader->draw(colorFrameBuffer->getTexture(), context->pixelRatio * rbbox, 1.0f, false);
		} else {
			matcapShader->draw(depthFrameBuffer->getTexture(), camera, context->pixelRatio * rbbox, context->getViewport(), surfaceColor.toRGBAf());
		}
		break;
	case 4:
		imageShader->draw(wireframeFrameBuffer->getTexture(), context->pixelRatio * rbbox, 1.0f, false);
		break;
	case 5:
		if (shadingType == 1) {
			imageShader->draw(colorFrameBuffer->getTexture(), context->pixelRatio * rbbox, 1.0f, false);
		} else {
			matcapShader->draw(depthFrameBuffer->getTexture(), camera, context->pixelRatio * rbbox, context->getViewport(), surfaceColor.toRGBAf());
		}
		break;
	default:
		break;
	}
	camera.setDirty(false);
	parametersDirty = false;
}

