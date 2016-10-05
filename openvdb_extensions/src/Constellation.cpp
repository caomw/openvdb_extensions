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
#include "Constellation.h"
#include "GeometryUtil.h"
#include "UtilitiesVDB.h"
#include <AlloyPLY.h>
#include <AlloyMath.h>
#include <AlloyMesh.h>
#include <list>
#include <openvdb/util/Util.h>
using namespace aly;
using namespace openvdb;
bool ReadConstellationFromFile(const std::string& file,Constellation& mesh){
	using namespace ply;
	int i, j, k;
	int numPts = 0, numPolys = 0;
	// open a PLY file for reading
	PLYReaderWriter ply;
	ply.openForReading(file);
	// Check to make sure that we can read geometry
	PlyElement *elem;
	int index;
	if ((elem = ply.findElement("vertex")) == nullptr
			||  ply.findProperty(elem, "x", &index) == nullptr
			||  ply.findProperty(elem, "y", &index) == nullptr
			||  ply.findProperty(elem, "z", &index) == nullptr
			|| (elem = ply.findElement("face")) == nullptr
			||  ply.findProperty(elem, "vertex_indices", &index) == nullptr) {
		std::cerr << "Cannot read geometry [" << file << "]" << std::endl;
		return false;
	}

	// Check for optional attribute data. We can handle intensity; and the
	// triplet red, green, blue.
	bool RGBPointsAvailable = false;
	bool hasParticleVelocity=false;
	bool hasNormals=false;
	bool hasVertexVelocity=false;
	//mesh.mTriIndexes.clear();
	//mesh.mQuadIndexes.clear();
	mesh.mQuads.clear();
	mesh.mTriangles.clear();
	mesh.mVertexes.clear();
	mesh.mParticles.clear();
	mesh.mParticleNormals.clear();
	mesh.mVertexNormals.clear();
	mesh.mColors.clear();
	mesh.mParticleVelocity.clear();
	mesh.mVertexVelocity.clear();
	if ((elem = ply.findElement("vertex")) != nullptr
			&&  ply.findProperty(elem, "red", &index) != nullptr
			&&  ply.findProperty(elem, "green", &index) != nullptr
			&&  ply.findProperty(elem, "blue", &index) != nullptr) {
		RGBPointsAvailable = true;
	}
	if ((elem = ply.findElement("vertex")) != nullptr
			&&  ply.findProperty(elem, "nx", &index) != nullptr
			&&  ply.findProperty(elem, "ny", &index) != nullptr
			&&  ply.findProperty(elem, "nz", &index) != nullptr) {
		hasNormals = true;
	}
	if ((elem = ply.findElement("vertex")) != nullptr
			&&  ply.findProperty(elem, "vx", &index) != nullptr
			&&  ply.findProperty(elem, "vy", &index) != nullptr
			&&  ply.findProperty(elem, "vz", &index) != nullptr) {
		hasVertexVelocity = true;
	}
	if ((elem = ply.findElement("face")) != nullptr&&  ply.findProperty(elem, "velocities", &index) != nullptr) {
		hasParticleVelocity = true;
	}
	int verts[256];
	float velocity[3];
	plyFace face;
	plyVertex vertex;
	face.verts=verts;
	face.velocity=velocity;
	memset(verts, 0, sizeof(verts));
	std::vector<std::string> elist = ply.getElementNames();
	std::string elemName;
	int numElems, nprops;
	// Okay, now we can grab the data
	for (i = 0; i < ply.getNumberOfElements(); i++) {
		//get the description of the first element */
		elemName = elist[i];
		ply.getElementDescription(elemName, &numElems, &nprops);
		// if we're on vertex elements, read them in
		if (elemName=="vertex") {
			// Create a list of points
			numPts = numElems;
			mesh.mVertexes.resize(numPts, Vec3s(0.0f));
			// Setup to read the PLY elements
			ply.getProperty( elemName, &MeshVertProps[0]);
			ply.getProperty( elemName, &MeshVertProps[1]);
			ply.getProperty( elemName, &MeshVertProps[2]);
			if (hasNormals) {
				mesh.mVertexNormals.resize(numPts);
				ply.getProperty( elemName, &MeshVertProps[3]);
				ply.getProperty( elemName, &MeshVertProps[4]);
				ply.getProperty( elemName, &MeshVertProps[5]);
			}
			if (hasVertexVelocity) {
				mesh.mVertexVelocity.resize(numPts);
				ply.getProperty( elemName, &MeshVertProps[6]);
				ply.getProperty( elemName, &MeshVertProps[7]);
				ply.getProperty( elemName, &MeshVertProps[8]);
			}
			if (RGBPointsAvailable) {
				mesh.mColors.resize(numPts);
				ply.getProperty( elemName, &MeshVertProps[9]);
				ply.getProperty( elemName, &MeshVertProps[10]);
				ply.getProperty( elemName, &MeshVertProps[11]);
			}
			for (j = 0; j < numPts; j++) {
				ply.getElement(&vertex);
				mesh.mVertexes[j] = Vec3s(vertex.x[0], vertex.x[1],
						vertex.x[2]);

				if (RGBPointsAvailable) {
					mesh.mColors[j] = Vec3s(vertex.red / 255.0f,
							vertex.green / 255.0f, vertex.blue / 255.0f);
				}
				if(hasNormals){
					mesh.mVertexNormals[j]=Vec3s(vertex.n[0],vertex.n[1],vertex.n[2]);
				}
				if(hasVertexVelocity){
					mesh.mVertexVelocity[j]=Vec3s(vertex.vel[0],vertex.vel[1],vertex.vel[2]);
				}
			}
		}			//if vertex
		else if (elemName=="face") {
			// Create a polygonal array
			numPolys = numElems;
			// Get the face properties
			ply.getProperty( elemName, &MeshFaceProps[0]);
			if(hasParticleVelocity){
				mesh.mParticleVelocity.resize(numPolys);
				ply.getProperty( elemName, &MeshFaceProps[1]);
			}
			for (j = 0; j < numPolys; j++) {
				ply.getElement(&face);
				if(hasParticleVelocity){
					Vec3s vel=Vec3s(face.velocity[0],face.velocity[1],face.velocity[2]);
					mesh.mParticleVelocity[j]=vel;
				}
				if (face.nverts == 4) {
					/*
					for (k = 0; k < face.nverts; k++) {
						mesh.mQuadIndexes.push_back(face.verts[k]);
					}
					*/
					mesh.mQuads.push_back(
							openvdb::Vec4I(face.verts[0], face.verts[1],
									face.verts[2], face.verts[3]));
				} else if (face.nverts == 3) {
					/*
					for (k = 0; k < face.nverts; k++) {
						mesh.mTriIndexes.push_back(face.verts[k]);
					}
					*/
					mesh.mTriangles.push_back(
							openvdb::Vec3I(face.verts[0], face.verts[1],
									face.verts[2]));
				}
			}
		}
	}
	if (mesh.mVertexes.size() > 0) {
		mesh.updateBoundingBox();
		return true;
	} else {
		return false;
	}
}
bool WriteConstellationToFile(const std::string& file,const Constellation& mesh) {
	using namespace ply;
	PLYReaderWriter ply;
	std::cout<<"Write "<<file<<std::endl;
	std::vector<std::string> elemNames = { "vertex", "face", "normal" };
	ply.openForWriting(file, elemNames, FileFormat::BINARY_LE );
	if (mesh.mVertexes.size() == 0)
		return false;
	int i, j, idx;
	bool usingTexture = (mesh.textureMap.size() > 0);
	// compute colors, if any
	int numPts = mesh.mVertexes.size();
	int numPolys = mesh.mQuads.size() + mesh.mTriangles.size();

	std::vector<unsigned char> pointColors;

	if (mesh.mColors.size() > 0) {
		size_t inc = 0;
		pointColors.resize(3 * mesh.mColors.size());
		for (i = 0; i < numPts; i++) {
			Vec3s d = mesh.mColors[i];
			pointColors[inc++] = (unsigned char) clamp(d[0] * 255.0f, 0.0f,
					255.0f);
			pointColors[inc++] = (unsigned char) clamp(d[1] * 255.0f, 0.0f,
					255.0f);
			pointColors[inc++] = (unsigned char) clamp(d[2] * 255.0f, 0.0f,
					255.0f);
		}
	}
	// describe what properties go into the vertex and face elements
	ply.elementCount("vertex", numPts);
	ply.describeProperty( "vertex", &MeshVertProps[0]);
	ply.describeProperty( "vertex", &MeshVertProps[1]);
	ply.describeProperty( "vertex", &MeshVertProps[2]);

	if(mesh.mVertexNormals.size()>0){
		ply.describeProperty( "vertex", &MeshVertProps[3]);
		ply.describeProperty( "vertex", &MeshVertProps[4]);
		ply.describeProperty( "vertex", &MeshVertProps[5]);
	}
	if (mesh.mVertexVelocity.size() > 0) {
		ply.describeProperty( "vertex", &MeshVertProps[6]);
		ply.describeProperty( "vertex", &MeshVertProps[7]);
		ply.describeProperty( "vertex", &MeshVertProps[8]);
	}
	if (mesh.mColors.size() > 0) {
		ply.describeProperty( "vertex", &MeshVertProps[9]);
		ply.describeProperty( "vertex", &MeshVertProps[10]);
		ply.describeProperty( "vertex", &MeshVertProps[11]);
	}
	ply.elementCount( "face", numPolys);

	if (usingTexture) {
		ply.describeProperty( "face", &MeshFaceProps[2]);
		ply.describeProperty( "face", &MeshFaceProps[3]);
		if(mesh.mParticleVelocity.size()>0)ply.describeProperty( "face", &MeshFaceProps[4]);
	} else {
		ply.describeProperty( "face", &MeshFaceProps[0]);
		if(mesh.mParticleVelocity.size()>0)ply.describeProperty( "face", &MeshFaceProps[1]);
	}

	// write a comment and an object information field
	ply.appendComment( "PLY File");
	if (usingTexture) {
		std::string comment = "TextureFile texture.png";
		ply.appendComment(comment);
	}
	ply.appendComment("ImageSci");

	// complete the header
	ply.headerComplete();
	// set up and write the vertex elements
	plyVertex vert;
	ply.putElementSetup("vertex");

	for (i = 0; i < numPts; i++) {
		Vec3s pt = mesh.mVertexes[i];
		vert.x[0] = pt[0];
		vert.x[1] = pt[1];
		vert.x[2] = pt[2];
		if(mesh.mVertexNormals.size()>0){
			Vec3s n=mesh.mVertexNormals[i];
			vert.n[0]=n[0];
			vert.n[1]=n[1];
			vert.n[2]=n[2];
		}
		if(mesh.mVertexVelocity.size()>0){
			Vec3s vel=mesh.mVertexVelocity[i];
			vert.vel[0]=vel[0];
			vert.vel[1]=vel[1];
			vert.vel[2]=vel[2];
		}
		if (pointColors.size() > 0) {
			idx = 3 * i;
			vert.red = pointColors[idx];
			vert.green = pointColors[idx + 1];
			vert.blue = pointColors[idx + 2];
		}
		ply.putElement((void *) &vert);
	}
	// set up and write the face elements
	plyFace face;
	plyFaceTexture faceT;
	int verts[256];
	Vec2s uvs[4];
	float vel[4];
	face.verts = verts;
	faceT.verts = verts;
	if(mesh.mParticleVelocity.size()>0){
		faceT.velocity=vel;
		face.velocity=vel;
		face.nvels=3;
		faceT.nvels=3;
	}
	faceT.uvs = (float*) uvs;
	ply.putElementSetup( "face");
	if (usingTexture) {
		int sz = mesh.mQuads.size();
		for (int i = 0; i < sz; i++) {
			faceT.nverts = 4;
			faceT.uvcount = 8;
			for (j = 0; j < 4; j++) {
				faceT.verts[j] = mesh.mQuads[i][j];
				uvs[j] = mesh.textureMap[4 * i + j];
			}
			if(faceT.velocity!=nullptr){
				Vec3s velocity=mesh.mParticleVelocity[i];
				faceT.velocity[0]=velocity[0];
				faceT.velocity[1]=velocity[1];
				faceT.velocity[2]=velocity[2];
			}
			ply.putElement((void *) &faceT);
		}
		sz = mesh.mTriangles.size();
		for (int i = 0; i < sz; i++) {
			faceT.nverts = 3;
			faceT.uvcount = 6;
			for (j = 0; j < 3; j++) {
				faceT.verts[j] = mesh.mTriangles[i][j];
				uvs[j] = mesh.textureMap[3 * i + j];
			}
			if(faceT.velocity!=nullptr){
				Vec3s velocity=mesh.mParticleVelocity[i];
				faceT.velocity[0]=velocity[0];
				faceT.velocity[1]=velocity[1];
				faceT.velocity[2]=velocity[2];
			}
			ply.putElement( (void *) &faceT);
		}
	} else {
		int sz = mesh.mQuads.size();
		for (int i = 0; i < sz; i++) {
			for (j = 0; j < 4; j++) {
				face.nverts = 4;
				face.verts[j] = mesh.mQuads[i][j];
			}
			if(faceT.velocity!=nullptr){
				Vec3s velocity=mesh.mParticleVelocity[i];
				faceT.velocity[0]=velocity[0];
				faceT.velocity[1]=velocity[1];
				faceT.velocity[2]=velocity[2];
			}
			ply.putElement( (void *) &face);
		}
		sz =mesh.mTriangles.size();
		for (int i = 0; i < sz; i++) {
			for (j = 0; j < 3; j++) {
				face.nverts = 3;
				face.verts[j] = mesh.mTriangles[i][j];
			}
			if(faceT.velocity!=nullptr){
				Vec3s velocity=mesh.mParticleVelocity[i];
				faceT.velocity[0]=velocity[0];
				faceT.velocity[1]=velocity[1];
				faceT.velocity[2]=velocity[2];
			}
			ply.putElement((void *) &face);
		}
	}
	return true;
}
Vec3s Constellation::closestPointOnEdge(const Vec3s& start,
		const SpringlNeighbor& ci) {
	Vec3s pt;
	Springl& nbr = springls[ci.springlId];
	DistanceToEdgeSqr(start, nbr[ci.edgeId], nbr[(ci.edgeId + 1) % nbr.size()],
			&pt);
	return pt;
}
void Constellation::copyTo(aly::Mesh& mesh) {
	Convert(mVertexes,mesh.vertexLocations);
	Convert(mQuads,mesh.quadIndexes);
	Convert(mTriangles,mesh.triIndexes);
	Convert(mVertexNormals,mesh.vertexNormals);
	mesh.updateBoundingBox();
	mesh.setDirty(true);
}
void Constellation::copyFrom(const Constellation& mesh) {
	size_t faceCount = mesh.mQuads.size()+mesh.mTriangles.size();
	size_t counter = 0;
	size_t pcounter = 0;
	springls.clear();
	mQuads.clear();
	mTriangles.clear();
	mVertexes.clear();
	mVertexes.resize(mesh.mQuads.size()*4 + mesh.mTriangles.size()*3);
	mParticles.resize(faceCount);
	mParticleNormals.resize(faceCount);
	mVertexNormals.resize(mVertexes.size());
	mParticleVelocity=mesh.mParticleVelocity;
	mVertexVelocity=mesh.mVertexVelocity;
	for (openvdb::Vec4I face : mesh.mQuads) {
		Springl springl(this);
		springl.offset = counter;
		springl.id = springls.size();
			mQuads.push_back(
					openvdb::Vec4I(counter, counter + 1, counter + 2,
							counter + 3));
			mVertexes[counter++] = mesh.mVertexes[face[0]];
			mVertexes[counter++] = mesh.mVertexes[face[1]];
			mVertexes[counter++] = mesh.mVertexes[face[2]];
			mVertexes[counter++] = mesh.mVertexes[face[3]];
			mParticles[pcounter] = springl.computeCentroid();
			openvdb::Vec3s norm = springl.computeNormal();
			mParticleNormals[pcounter] = norm;
			mVertexNormals[counter - 1] = norm;
			mVertexNormals[counter - 2] = norm;
			mVertexNormals[counter - 3] = norm;
			mVertexNormals[counter - 4] = norm;
			springls.push_back(springl);

		pcounter++;
	}
	for (openvdb::Vec3I face : mesh.mTriangles) {
			Springl springl(this);
			springl.offset = counter;
			springl.id = springls.size();
				mTriangles.push_back(openvdb::Vec3I(counter, counter + 1, counter + 2));
				mVertexes[counter++] = mesh.mVertexes[face[0]];
				mVertexes[counter++] = mesh.mVertexes[face[1]];
				mVertexes[counter++] = mesh.mVertexes[face[2]];
				mParticles[pcounter] = springl.computeCentroid();
				openvdb::Vec3s norm = springl.computeNormal();
				mParticleNormals[pcounter] = norm;
				mVertexNormals[counter - 1] = norm;
				mVertexNormals[counter - 2] = norm;
				mVertexNormals[counter - 3] = norm;
				springls.push_back(springl);
			pcounter++;
		}
	updateBoundingBox();
}
Constellation::Constellation() : mPose(openvdb::math::Mat4f::identity()){
}
void Constellation::reset() {
	mVertexes.clear();
	mVertexNormals.clear();
	mParticleNormals.clear();
	mParticles.clear();
	mQuads.clear();
	mTriangles.clear();
	mVertexAuxBuffer.clear();
}

void Constellation::initialize() {
	size_t faceCount = mQuads.size()+mTriangles.size();
	size_t counter = 0;
	size_t pcounter = 0;
	springls.clear();
	mQuads.clear();
	mTriangles.clear();
	mVertexes.clear();
	mVertexes.resize(mQuads.size()*4 + mTriangles.size()*3);
	mParticles.resize(faceCount);
	mParticleNormals.resize(faceCount);
	mVertexNormals.resize(mVertexes.size());
	for (openvdb::Vec3I face : mTriangles) {
		Springl springl(this);
		springl.offset = counter;
		springl.id = springls.size();
			mTriangles.push_back(openvdb::Vec3I(counter, counter + 1, counter + 2));
			mVertexes[counter++] = mVertexes[face[0]];
			mVertexes[counter++] = mVertexes[face[1]];
			mVertexes[counter++] = mVertexes[face[2]];
			mParticles[pcounter] = springl.computeCentroid();
			openvdb::Vec3s norm = springl.computeNormal();
			mParticleNormals[pcounter] = norm;
			mVertexNormals[counter - 1] = norm;
			mVertexNormals[counter - 2] = norm;
			mVertexNormals[counter - 3] = norm;

			springls.push_back(springl);

		pcounter++;
	}
	for (openvdb::Vec4I face : mQuads) {
			Springl springl(this);
			springl.offset = counter;
			springl.id = springls.size();
				mQuads.push_back(
						openvdb::Vec4I(counter, counter + 1, counter + 2,
								counter + 3));
				mVertexes[counter++] = mVertexes[face[0]];
				mVertexes[counter++] = mVertexes[face[1]];
				mVertexes[counter++] = mVertexes[face[2]];
				mVertexes[counter++] = mVertexes[face[3]];
				mParticles[pcounter] = springl.computeCentroid();
				openvdb::Vec3s norm = springl.computeNormal();
				mParticleNormals[pcounter] = norm;
				mVertexNormals[counter - 1] = norm;
				mVertexNormals[counter - 2] = norm;
				mVertexNormals[counter - 3] = norm;
				mVertexNormals[counter - 4] = norm;
				springls.push_back(springl);
			pcounter++;
		}
	updateBoundingBox();
}
void Constellation::dilate(float distance) {
	int vertCount = mVertexes.size();
	Vec3s norm;
	for (int i = 0; i < vertCount; i++) {
		norm = mVertexNormals[i];
		mVertexes[i] += norm * distance;
	}
}
void Constellation::updateVertexNormals(int SMOOTH_ITERATIONS, float DOT_TOLERANCE) {

	Index32 sz = mTriangles.size();
	Vec3s pt;
	mVertexNormals.resize(mVertexes.size(), Vec3f(0.0f));
	for (Index32 i = 0; i < sz; i ++) {
		Vec3s v1 = mVertexes[mTriangles[i][0]];
		Vec3s v2 = mVertexes[mTriangles[i][1]];
		Vec3s v3 = mVertexes[mTriangles[i][2]];
		Vec3f norm = (v2 - v1).cross(v3 - v1);
		mVertexNormals[mTriangles[i][0]] += norm;
		mVertexNormals[mTriangles[i][1]] += norm;
		mVertexNormals[mTriangles[i][2]] += norm;
	}
	sz = mQuads.size();
	for (int i = 0; i < (int)sz; i ++) {
		Vec3s v1 = mVertexes[mQuads[i][0]];
		Vec3s v2 = mVertexes[mQuads[i][1]];
		Vec3s v3 = mVertexes[mQuads[i][2]];
		Vec3s v4 = mVertexes[mQuads[i][3]];
		Vec3f norm = (v1 - pt).cross(v2 - pt);
		norm += (v2 - pt).cross(v3 - pt);
		norm += (v3 - pt).cross(v4 - pt);
		norm += (v4 - pt).cross(v1 - pt);
		mVertexNormals[mQuads[i][0]] += norm;
		mVertexNormals[mQuads[i][1]] += norm;
		mVertexNormals[mQuads[i][2]] += norm;
		mVertexNormals[mQuads[i][3]] += norm;
	}
#pragma omp parallel for
	for (size_t n=0;n<mVertexNormals.size();n++) {
		mVertexNormals[n].normalize(1E-6f);
	}
	if (SMOOTH_ITERATIONS > 0) {
		int vertCount = mVertexes.size();
		std::vector<Vec3f> tmp(vertCount);
		std::vector<std::list<int>> vertNbrs(vertCount);
		int indexCount = mQuads.size();
		for (int i = 0; i < indexCount; i ++) {
			int v1 =mQuads[i][0];
			int v2 =mQuads[i][1];
			int v3 =mQuads[i][2];
			int v4 =mQuads[i][3];

			vertNbrs[v1].push_back(v2);
			vertNbrs[v2].push_back(v3);
			vertNbrs[v3].push_back(v1);

			vertNbrs[v3].push_back(v4);
			vertNbrs[v4].push_back(v1);
			vertNbrs[v1].push_back(v3);
		}
		for (int iter = 0; iter < SMOOTH_ITERATIONS; iter++) {
#pragma omp parallel for
			for (int i = 0; i < vertCount; i++) {
				Vec3f norm = mVertexNormals[i];
				Vec3f avg = Vec3f(0.0f);
				for (int nbr : vertNbrs[i]) {
					Vec3s nnorm = mVertexNormals[nbr];
					;
					if (norm.dot(nnorm) > DOT_TOLERANCE) {
						avg += nnorm;
					} else {
						avg += norm;
					}
				}
				avg.normalize();
				tmp[i] = avg;
			}
			mVertexNormals = tmp;
		}
	}
}
float Constellation::estimateVoxelSize(int stride) {
	int count = 0;
	//float maxLength = 0.0f;
	int sz = mTriangles.size();
	float mEstimatedVoxelSize = 0.0f;
	for (int i = 0; i < sz; i += stride) {
		Vec3s v1 = mVertexes[mTriangles[i][0]];
		Vec3s v2 = mVertexes[mTriangles[i][1]];
		Vec3s v3 = mVertexes[mTriangles[i][2]];
		float e1 = (v1 - v2).length();
		float e2 = (v1 - v3).length();
		float e3 = (v2 - v3).length();
		//maxLength = std::max(std::max(e1, e2), std::max(maxLength, e3));
		mEstimatedVoxelSize += e1 + e2 + e3;
	}
	count = sz / stride;
	sz = mQuads.size();
	for (int i = 0; i < sz; i += stride) {
		Vec3s v1 = mVertexes[mQuads[i][0]];
		Vec3s v2 = mVertexes[mQuads[i][1]];
		Vec3s v3 = mVertexes[mQuads[i][2]];
		Vec3s v4 = mVertexes[mQuads[i][3]];
		float e1 = (v1 - v2).length();
		float e2 = (v2 - v3).length();
		float e3 = (v3 - v4).length();
		float e4 = (v4 - v1).length();
		//maxLength = std::max(maxLength,std::max(std::max(e1, e2), std::max(e3, e4)));
		mEstimatedVoxelSize += e1 + e2 + e3 + e4;
	}
	count += sz / stride;
	mEstimatedVoxelSize /= count;

	std::cout << "Estimated voxel size =" << mEstimatedVoxelSize << std::endl;
	return mEstimatedVoxelSize;
}
openvdb::math::BBox<openvdb::Vec3d>& Constellation::updateBoundingBox() {
	const int BATCHES = 32;
	Vec3s minPt(std::numeric_limits<float>::max(),
			std::numeric_limits<float>::max(),
			std::numeric_limits<float>::max());
	std::vector<Vec3s> minPtBatch(BATCHES,
			Vec3s(std::numeric_limits<float>::max(),
					std::numeric_limits<float>::max(),
					std::numeric_limits<float>::max()));
	Vec3s maxPt(std::numeric_limits<float>::min(),
			std::numeric_limits<float>::min(),
			std::numeric_limits<float>::min());
	std::vector<Vec3s> maxPtBatch(BATCHES,
			Vec3s(std::numeric_limits<float>::min(),
					std::numeric_limits<float>::min(),
					std::numeric_limits<float>::min()));
	int SZ = mVertexes.size();
	int batchSize = (SZ % BATCHES == 0) ? SZ / BATCHES : SZ / BATCHES + 1;
#pragma omp for
	for (int b = 0; b < BATCHES; b++) {
		int sz = std::min(SZ, batchSize * (b + 1));
		for (int idx = b * batchSize; idx < (int)sz; idx++) {
			Vec3s& pt = mVertexes[idx];
			minPtBatch[b][0] = std::min(minPtBatch[b][0], pt[0]);
			minPtBatch[b][1] = std::min(minPtBatch[b][1], pt[1]);
			minPtBatch[b][2] = std::min(minPtBatch[b][2], pt[2]);

			maxPtBatch[b][0] = std::max(maxPtBatch[b][0], pt[0]);
			maxPtBatch[b][1] = std::max(maxPtBatch[b][1], pt[1]);
			maxPtBatch[b][2] = std::max(maxPtBatch[b][2], pt[2]);
		}
	}

	for (int b = 0; b < BATCHES; b++) {
		minPt[0] = std::min(minPtBatch[b][0], minPt[0]);
		minPt[1] = std::min(minPtBatch[b][1], minPt[1]);
		minPt[2] = std::min(minPtBatch[b][2], minPt[2]);

		maxPt[0] = std::max(maxPtBatch[b][0], maxPt[0]);
		maxPt[1] = std::max(maxPtBatch[b][1], maxPt[1]);
		maxPt[2] = std::max(maxPtBatch[b][2], maxPt[2]);
	}
	mBoundingBox = openvdb::math::BBox<openvdb::Vec3d>(minPt, maxPt);
	return mBoundingBox;
}
void Constellation::scale(float sc) {
#pragma omp for
	for (size_t i = 0; i < mVertexes.size(); i++) {
		mVertexes[i] *= sc;
	}
	mBoundingBox.min() *= static_cast<double>(sc);
	mBoundingBox.max() *= static_cast<double>(sc);

}
void Constellation::mapIntoBoundingBox(float voxelSize) {
	Vec3s minPt = mBoundingBox.min();
#pragma omp for
	for (size_t i = 0; i < mVertexes.size(); i++) {
		Vec3s& pt = mVertexes[i];
		pt = (pt - minPt) / voxelSize;
	}
}
void Constellation::mapOutOfBoundingBox(float voxelSize) {
	Vec3s minPt = mBoundingBox.min();
#pragma omp for
	for (size_t i = 0; i < mVertexes.size(); i++) {
		Vec3s& pt = mVertexes[i];
		pt = pt * voxelSize + minPt;
	}
}

Constellation::~Constellation() {
	// TODO Auto-generated destructor stub
}




