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
#include "UtilitiesVDB.h"
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/Prune.h>
#include <openvdb/tree/LeafManager.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/Dense.h>
#include <AlloyMeshPrimitives.h>
#include <AlloyPLY.h>
using namespace openvdb;
using namespace aly;

bool WriteToRawFile(const std::string& file,const RegularGrid<float>& dense) {
	std::ostringstream vstr;
	std::string fileName=GetFileWithoutExtension(file);
	vstr << fileName << ".raw";

	FILE* f = fopen(vstr.str().c_str(), "wb");
	openvdb::CoordBBox bbox = dense.bbox();
	std::cout << "Grid size " << dense.valueCount() << std::endl;
	openvdb::Coord dims = bbox.max() - bbox.min() + openvdb::Coord(1, 1, 1);
	std::cout << "Dimensions " << dims << std::endl;
	openvdb::Coord P(0, 0, 0);
	for (P[2] = bbox.min()[2]; P[2] <= bbox.max()[2]; ++P[2]) {
		for (P[1] = bbox.min()[1]; P[1] <= bbox.max()[1]; ++P[1]) {

			for (P[0] = bbox.min()[0]; P[0] <= bbox.max()[0]; ++P[0]) {
				float val = dense.getValue(P);
				fwrite(&val, sizeof(float), 1, f);
			}
		}
	}
	fclose(f);
	std::stringstream xmlFile;
	xmlFile << fileName << ".xml";
	std::cout <<"Saving "<< xmlFile.str() <<" ...";
	std::stringstream sstr;
	sstr << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	sstr << "<!-- MIPAV header file -->\n";
	sstr<< "<image xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" nDimensions=\"3\">\n";
	sstr << "	<Dataset-attributes>\n";
	sstr << "		<Image-offset>0</Image-offset>\n";
	sstr << "		<Data-type>Float</Data-type>\n";
	sstr << "		<Endianess>Little</Endianess>\n";
	sstr << "		<Extents>" << dims[0] << "</Extents>\n";
	sstr << "		<Extents>" << dims[1] << "</Extents>\n";
	sstr << "		<Extents>" << dims[2] << "</Extents>\n";
	sstr << "		<Resolutions>\n";
	sstr << "			<Resolution>1.0</Resolution>\n";
	sstr << "			<Resolution>1.0</Resolution>\n";
	sstr << "			<Resolution>1.0</Resolution>\n";
	sstr << "		</Resolutions>\n";
	sstr << "		<Slice-spacing>1.0</Slice-spacing>\n";
	sstr << "		<Slice-thickness>0.0</Slice-thickness>\n";
	sstr << "		<Units>Millimeters</Units>\n";
	sstr << "		<Units>Millimeters</Units>\n";
	sstr << "		<Units>Millimeters</Units>\n";
	sstr << "		<Compression>none</Compression>\n";
	sstr << "		<Orientation>Unknown</Orientation>\n";
	sstr << "		<Subject-axis-orientation>Unknown</Subject-axis-orientation>\n";
	sstr << "		<Subject-axis-orientation>Unknown</Subject-axis-orientation>\n";
	sstr << "		<Subject-axis-orientation>Unknown</Subject-axis-orientation>\n";
	sstr << "		<Origin>0.0</Origin>\n";
	sstr << "		<Origin>0.0</Origin>\n";
	sstr << "		<Origin>0.0</Origin>\n";
	sstr << "		<Modality>Unknown Modality</Modality>\n";
	sstr << "	</Dataset-attributes>\n";
	sstr << "</image>\n";
	std::ofstream myfile;
	myfile.open(xmlFile.str().c_str(), std::ios_base::out);
	myfile << sstr.str();
	myfile.close();
	std::cout<<" done."<<std::endl;
	return true;
}
void ConvertLevelSetToConstellation(const openvdb::FloatGrid::Ptr& grid,openvdb::tools::VolumeToMesh& mesher,Constellation& mesh) {
	// Copy points and generate point normals.
	openvdb::math::GenericMap map(grid->transform());
	mesh.mVertexes.clear();
	mesh.mVertexNormals.clear();
	mesh.mQuads.clear();
	mesh.mTriangles.clear();
	mesher(*grid);
	mesh.mVertexes.resize(mesher.pointListSize());
	Index64 N = mesher.pointListSize();
	for (Index64 n = 0; n < N; ++n) {
		mesh.mVertexes[n] = mesher.pointList()[n];	//map.applyInverseMap(
	}
	openvdb::tools::PolygonPoolList& polygonPoolList = mesher.polygonPoolList();
	for (Index64 n = 0, N = mesher.polygonPoolListSize(); n < N; ++n) {
		const openvdb::tools::PolygonPool& polygons = polygonPoolList[n];
		for (Index64 i = 0, I = polygons.numQuads(); i < I; ++i) {
			const openvdb::Vec4I& quad = polygons.quad(i);
			mesh.mQuads.push_back(
					openvdb::Vec4I(quad[3], quad[2], quad[1], quad[0]));

		}
		for (Index64 i = 0, I = polygons.numTriangles(); i < I; ++i) {
			const openvdb::Vec3I& quad = polygons.triangle(i);
			mesh.mTriangles.push_back(
					openvdb::Vec3I(quad[2], quad[1], quad[0]));
		}
	}
	mesh.updateVertexNormals(4);
	mesh.updateBoundingBox();
}
void ConvertLevelSetToConstellation(const FloatGrid::Ptr& grid,Constellation& mesh) {
	openvdb::tools::VolumeToMesh mesher(0.0);
	ConvertLevelSetToConstellation(grid,mesher,mesh);
}
bool WriteToRawFile(const std::string& fileName,const openvdb::FloatGrid::Ptr& grid) {
	using namespace openvdb::tools;
	std::ostringstream vstr;
	std::string file=GetFileWithoutExtension(fileName);
	vstr << file << ".raw";
	FILE* f = fopen(vstr.str().c_str(), "wb");
	openvdb::CoordBBox bbox = grid->evalActiveVoxelBoundingBox();
	Dense<float> dense(bbox); //LayoutZYX is the default
	copyToDense(*grid, dense);
	Coord dims = bbox.max() - bbox.min() + Coord(1, 1, 1);
	openvdb::Coord P(0, 0, 0);
	for (P[2] = bbox.min()[2]; P[2] <= bbox.max()[2]; ++P[2]) {
		for (P[1] = bbox.min()[1]; P[1] <= bbox.max()[1]; ++P[1]) {
			for (P[0] = bbox.min()[0]; P[0] <= bbox.max()[0]; ++P[0]) {
				float val = dense.getValue(P);
				fwrite(&val, sizeof(float), 1, f);
			}
		}
	}
	fclose(f);
	std::stringstream sstr;
	sstr << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	sstr << "<!-- MIPAV header file -->\n";
	sstr
			<< "<image xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" nDimensions=\"3\">\n";
	sstr << "	<Dataset-attributes>\n";
	sstr << "		<Image-offset>0</Image-offset>\n";
	sstr << "		<Data-type>Float</Data-type>\n";
	sstr << "		<Endianess>Little</Endianess>\n";
	sstr << "		<Extents>" << dims[0] << "</Extents>\n";
	sstr << "		<Extents>" << dims[1] << "</Extents>\n";
	sstr << "		<Extents>" << dims[2] << "</Extents>\n";
	sstr << "		<Resolutions>\n";
	sstr << "			<Resolution>1.0</Resolution>\n";
	sstr << "			<Resolution>1.0</Resolution>\n";
	sstr << "			<Resolution>1.0</Resolution>\n";
	sstr << "		</Resolutions>\n";
	sstr << "		<Slice-spacing>1.0</Slice-spacing>\n";
	sstr << "		<Slice-thickness>0.0</Slice-thickness>\n";
	sstr << "		<Units>Millimeters</Units>\n";
	sstr << "		<Units>Millimeters</Units>\n";
	sstr << "		<Units>Millimeters</Units>\n";
	sstr << "		<Compression>none</Compression>\n";
	sstr << "		<Orientation>Unknown</Orientation>\n";
	sstr << "		<Subject-axis-orientation>Unknown</Subject-axis-orientation>\n";
	sstr << "		<Subject-axis-orientation>Unknown</Subject-axis-orientation>\n";
	sstr << "		<Subject-axis-orientation>Unknown</Subject-axis-orientation>\n";
	sstr << "		<Origin>0.0</Origin>\n";
	sstr << "		<Origin>0.0</Origin>\n";
	sstr << "		<Origin>0.0</Origin>\n";
	sstr << "		<Modality>Unknown Modality</Modality>\n";
	sstr << "	</Dataset-attributes>\n";
	sstr << "</image>\n";
	std::ofstream myfile;
	std::stringstream xmlFile;
	xmlFile <<file << ".xml";
	myfile.open(xmlFile.str().c_str(), std::ios_base::out);
	myfile << sstr.str();
	myfile.close();
	std::cout << xmlFile.str() << std::endl;
	return true;
}
bool WriteToRawFile( const std::string& file,const MACGrid<float>& mac) {
	std::string fileName=GetFileWithoutExtension(file);
	bool r1 = WriteToRawFile(fileName + "_x.xml",mac[0]);
	bool r2 = WriteToRawFile(fileName + "_y.xml",mac[1]);
	bool r3 = WriteToRawFile(fileName + "_z.xml",mac[2]);
	return (r1 && r2 && r3);
}
bool WriteToRawFile(const std::string& fileName,const openvdb::VectorGrid::Ptr& grid) {
	using namespace openvdb::tools;
	std::ostringstream vstr;
	vstr << fileName << ".raw";
	FILE* f = fopen(vstr.str().c_str(), "wb");
	openvdb::CoordBBox bbox = grid->evalActiveVoxelBoundingBox();
	Dense<Vec3f> dense(bbox); //LayoutZYX is the default
	copyToDense(*grid, dense);
	Coord dims = bbox.max() - bbox.min() + Coord(1, 1, 1);
	openvdb::Coord P(0, 0, 0);
	for (P[2] = bbox.min()[2]; P[2] <= bbox.max()[2]; ++P[2]) {
		for (P[1] = bbox.min()[1]; P[1] <= bbox.max()[1]; ++P[1]) {
			for (P[0] = bbox.min()[0]; P[0] <= bbox.max()[0]; ++P[0]) {
				Vec3f val = dense.getValue(P);
				fwrite(&val[0], sizeof(float), 1, f);
			}
		}
	}

	for (P[2] = bbox.min()[2]; P[2] <= bbox.max()[2]; ++P[2]) {
		for (P[1] = bbox.min()[1]; P[1] <= bbox.max()[1]; ++P[1]) {
			for (P[0] = bbox.min()[0]; P[0] <= bbox.max()[0]; ++P[0]) {
				Vec3f val = dense.getValue(P);
				fwrite(&val[1], sizeof(float), 1, f);
			}
		}
	}

	for (P[2] = bbox.min()[2]; P[2] <= bbox.max()[2]; ++P[2]) {
		for (P[1] = bbox.min()[1]; P[1] <= bbox.max()[1]; ++P[1]) {
			for (P[0] = bbox.min()[0]; P[0] <= bbox.max()[0]; ++P[0]) {
				Vec3f val = dense.getValue(P);
				fwrite(&val[2], sizeof(float), 1, f);
			}
		}
	}

	fclose(f);
	std::stringstream sstr;
	sstr << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	sstr << "<!-- MIPAV header file -->\n";
	sstr
			<< "<image xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" nDimensions=\"4\">\n";
	sstr << "	<Dataset-attributes>\n";
	sstr << "		<Image-offset>0</Image-offset>\n";
	sstr << "		<Data-type>Float</Data-type>\n";
	sstr << "		<Endianess>Little</Endianess>\n";
	sstr << "		<Extents>" << dims[0] << "</Extents>\n";
	sstr << "		<Extents>" << dims[1] << "</Extents>\n";
	sstr << "		<Extents>" << dims[2] << "</Extents>\n";
	sstr << "		<Extents>3</Extents>\n";
	sstr << "		<Resolutions>\n";
	sstr << "			<Resolution>1.0</Resolution>\n";
	sstr << "			<Resolution>1.0</Resolution>\n";
	sstr << "			<Resolution>1.0</Resolution>\n";
	sstr << "		</Resolutions>\n";
	sstr << "		<Slice-spacing>1.0</Slice-spacing>\n";
	sstr << "		<Slice-thickness>0.0</Slice-thickness>\n";
	sstr << "		<Units>Millimeters</Units>\n";
	sstr << "		<Units>Millimeters</Units>\n";
	sstr << "		<Units>Millimeters</Units>\n";
	sstr << "		<Compression>none</Compression>\n";
	sstr << "		<Orientation>Unknown</Orientation>\n";
	sstr << "		<Subject-axis-orientation>Unknown</Subject-axis-orientation>\n";
	sstr << "		<Subject-axis-orientation>Unknown</Subject-axis-orientation>\n";
	sstr << "		<Subject-axis-orientation>Unknown</Subject-axis-orientation>\n";
	sstr << "		<Origin>0.0</Origin>\n";
	sstr << "		<Origin>0.0</Origin>\n";
	sstr << "		<Origin>0.0</Origin>\n";
	sstr << "		<Modality>Unknown Modality</Modality>\n";
	sstr << "	</Dataset-attributes>\n";
	sstr << "</image>\n";
	std::ofstream myfile;
	std::stringstream xmlFile;
	xmlFile << fileName << ".xml";
	myfile.open(xmlFile.str().c_str(), std::ios_base::out);
	myfile << sstr.str();
	myfile.close();
	std::cout << xmlFile.str() << std::endl;
	return true;
}
bool WriteToRawFile( const std::string& fileName,const openvdb::Int32Grid::Ptr& grid) {
	std::ostringstream vstr;
	using namespace openvdb::tools;
	std::string file=GetFileWithoutExtension(fileName);
	vstr << file << ".raw";
	FILE* f = fopen(vstr.str().c_str(), "wb");
	openvdb::CoordBBox bbox = grid->evalActiveVoxelBoundingBox();
	Dense<Index32> dense(bbox); //LayoutZYX is the default
	copyToDense(*grid, dense);
	Coord dims = bbox.max() - bbox.min() + Coord(1, 1, 1);
	openvdb::Coord P(0, 0, 0);
	for (P[2] = bbox.min()[2]; P[2] <= bbox.max()[2]; ++P[2]) {
		for (P[1] = bbox.min()[1]; P[1] <= bbox.max()[1]; ++P[1]) {

			for (P[0] = bbox.min()[0]; P[0] <= bbox.max()[0]; ++P[0]) {
				Index32 val = dense.getValue(P);
				fwrite(&val, sizeof(Index32), 1, f);
			}
		}
	}
	fclose(f);
	std::stringstream sstr;
	sstr << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	sstr << "<!-- MIPAV header file -->\n";
	sstr
			<< "<image xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" nDimensions=\"3\">\n";
	sstr << "	<Dataset-attributes>\n";
	sstr << "		<Image-offset>0</Image-offset>\n";
	sstr << "		<Data-type>Unsigned Integer</Data-type>\n";
	sstr << "		<Endianess>Little</Endianess>\n";
	sstr << "		<Extents>" << dims[0] << "</Extents>\n";
	sstr << "		<Extents>" << dims[1] << "</Extents>\n";
	sstr << "		<Extents>" << dims[2] << "</Extents>\n";
	sstr << "		<Resolutions>\n";
	sstr << "			<Resolution>1.0</Resolution>\n";
	sstr << "			<Resolution>1.0</Resolution>\n";
	sstr << "			<Resolution>1.0</Resolution>\n";
	sstr << "		</Resolutions>\n";
	sstr << "		<Slice-spacing>1.0</Slice-spacing>\n";
	sstr << "		<Slice-thickness>0.0</Slice-thickness>\n";
	sstr << "		<Units>Millimeters</Units>\n";
	sstr << "		<Units>Millimeters</Units>\n";
	sstr << "		<Units>Millimeters</Units>\n";
	sstr << "		<Compression>none</Compression>\n";
	sstr << "		<Orientation>Unknown</Orientation>\n";
	sstr << "		<Subject-axis-orientation>Unknown</Subject-axis-orientation>\n";
	sstr << "		<Subject-axis-orientation>Unknown</Subject-axis-orientation>\n";
	sstr << "		<Subject-axis-orientation>Unknown</Subject-axis-orientation>\n";
	sstr << "		<Origin>0.0</Origin>\n";
	sstr << "		<Origin>0.0</Origin>\n";
	sstr << "		<Origin>0.0</Origin>\n";
	sstr << "		<Modality>Unknown Modality</Modality>\n";
	sstr << "	</Dataset-attributes>\n";
	sstr << "</image>\n";
	std::ofstream myfile;
	std::stringstream xmlFile;
	xmlFile << file << ".xml";
	myfile.open(xmlFile.str().c_str(), std::ios_base::out);
	myfile << sstr.str();
	myfile.close();
	std::cout << xmlFile.str() << std::endl;
	return true;
}
void ConvertLevelSetToVolume(const openvdb::FloatGrid::Ptr& grid, aly::Volume1f& volume) {
	openvdb::CoordBBox bbox = grid->evalActiveVoxelBoundingBox();
	openvdb::tools::Dense<float> dense(bbox);
	openvdb::tools::copyToDense(*grid, dense);
	Coord dims = bbox.max() - bbox.min() + Coord(1, 1, 1);
	openvdb::Coord P(0, 0, 0);
	size_t index=0;
	volume.resize(dims[0],dims[1],dims[2]);
	for (P[2] = bbox.min()[2]; P[2] <= bbox.max()[2]; ++P[2]) {
		for (P[1] = bbox.min()[1]; P[1] <= bbox.max()[1]; ++P[1]) {
			for (P[0] = bbox.min()[0]; P[0] <= bbox.max()[0]; ++P[0]) {
				volume[index++] = float1(dense.getValue(P));
			}
		}
	}
	Coord minPt=bbox.min();
	volume.setPosition(int3(minPt[0],minPt[1],minPt[2]));
}
void ConvertVectorFieldToVolume(const openvdb::VectorGrid::Ptr& grid, aly::Volume3f& volume) {
	openvdb::CoordBBox bbox = grid->evalActiveVoxelBoundingBox();
	openvdb::tools::Dense<Vec3f> dense(bbox);
	openvdb::tools::copyToDense(*grid, dense);
	Coord dims = bbox.max() - bbox.min() + Coord(1, 1, 1);
	openvdb::Coord P(0, 0, 0);
	size_t index=0;
	volume.resize(dims[0],dims[1],dims[2]);
	for (P[2] = bbox.min()[2]; P[2] <= bbox.max()[2]; ++P[2]) {
		for (P[1] = bbox.min()[1]; P[1] <= bbox.max()[1]; ++P[1]) {
			for (P[0] = bbox.min()[0]; P[0] <= bbox.max()[0]; ++P[0]) {
				Vec3f val = dense.getValue(P);
				volume[index++] = float3(val[0],val[1],val[2]);
			}
		}
	}
	Coord minPt=bbox.min();
	volume.setPosition(int3(minPt[0],minPt[1],minPt[2]));
}
void ConvertLevelSetToBBoxTree(const openvdb::FloatGrid::Ptr& grid, aly::Mesh& mesh, int minDepth, int maxDepth) {
	using openvdb::Index64;
	openvdb::Vec3d ptn;
	float4 color;
	openvdb::CoordBBox bbox;
	Index64 idx = 0;
	mesh.clear();
	for (auto iter = grid->tree().cbeginNode(); iter; ++iter) {
		const int level = iter.getLevel();
		if (level >= minDepth && level <= maxDepth) {
			iter.getBoundingBox(bbox);
			float3 minPt(bbox.min().x() - 0.5, bbox.min().y() - 0.5, bbox.min().z() - 0.5);
			float3 maxPt(bbox.max().x() + 0.5, bbox.max().y() + 0.5, bbox.max().z() + 0.5);
			mesh.vertexLocations.push_back(float3(minPt[0], minPt[1], minPt[2]));
			mesh.vertexLocations.push_back(float3(maxPt[0], minPt[1], minPt[2]));
			mesh.vertexLocations.push_back(float3(maxPt[0], maxPt[1], minPt[2]));
			mesh.vertexLocations.push_back(float3(minPt[0], maxPt[1], minPt[2]));
			mesh.vertexLocations.push_back(float3(minPt[0], minPt[1], maxPt[2]));
			mesh.vertexLocations.push_back(float3(maxPt[0], minPt[1], maxPt[2]));
			mesh.vertexLocations.push_back(float3(maxPt[0], maxPt[1], maxPt[2]));
			mesh.vertexLocations.push_back(float3(minPt[0], maxPt[1], maxPt[2]));
			uint4 off(idx);
			mesh.quadIndexes.push_back(uint4(3, 2, 1, 0) + off);
			mesh.quadIndexes.push_back(uint4(6, 7, 4, 5) + off);
			mesh.quadIndexes.push_back(uint4(5, 1, 2, 6) + off);
			mesh.quadIndexes.push_back(uint4(0, 4, 7, 3) + off);
			mesh.quadIndexes.push_back(uint4(6, 2, 3, 7) + off);
			mesh.quadIndexes.push_back(uint4(4, 0, 1, 5) + off);
			for (Index64 n = 0; n < 8; ++n) {
				mesh.vertexColors.push_back(color);
			}
			idx += 8;
		}
	}
	mesh.updateVertexNormals();
	mesh.updateBoundingBox();
	mesh.setDirty(true);
}

void ConvertLevelSetToMesh(const aly::Mesh& mesh, openvdb::FloatGrid::Ptr& grid) {
	openvdb::math::Transform::Ptr trans = openvdb::math::Transform::createLinearTransform(1.0);
	std::vector<openvdb::Vec3s> vertexes(mesh.vertexLocations.size());
	std::vector<openvdb::Vec4I> quads;
	std::vector<openvdb::Vec3I> tris;
	for (uint4 quad : mesh.quadIndexes.data) {
		quads.push_back(Vec4I(quad.x, quad.y, quad.z, quad.w));
	}
	for (uint3 tri : mesh.triIndexes.data) {
		tris.push_back(Vec3I(tri.x, tri.y, tri.z));
	}
	std::memcpy(vertexes.data(), mesh.vertexLocations.ptr(), mesh.vertexLocations.size() * sizeof(Vec3s));
	grid=openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(*trans,vertexes,tris,quads, float(LEVEL_SET_HALF_WIDTH));

}

void ConvertLevelSetToMesh(const openvdb::FloatGrid::Ptr& grid, aly::Mesh& mesh) {
	openvdb::tools::VolumeToMesh mesher(0.0);
	mesher(*grid);
	openvdb::math::GenericMap map(grid->transform());
	mesh.vertexLocations.clear();
	mesh.vertexNormals.clear();
	mesh.triIndexes.clear();
	mesh.quadIndexes.clear();
	Index64 N = mesher.pointListSize();
	mesh.vertexLocations.resize(N);
	boost::scoped_array<openvdb::Vec3s>& points = mesher.pointList();
	for (Index64 n = 0; n < N; ++n) {
		Vec3d pt = map.applyInverseMap(points[n]);
		mesh.vertexLocations[n] = float3(pt[0], pt[1], pt[2]);
	}
	openvdb::tools::PolygonPoolList& polygonPoolList = mesher.polygonPoolList();
	N = mesher.polygonPoolListSize();
	for (Index64 n = 0; n < N; ++n) {
		const openvdb::tools::PolygonPool& polygons = polygonPoolList[n];
		for (Index64 i = 0, I = polygons.numQuads(); i < I; ++i) {
			const openvdb::Vec4I& quad = polygons.quad(i);
			mesh.quadIndexes.push_back(uint4(quad[3], quad[2], quad[1], quad[0]));
		}
		for (Index64 i = 0, I = polygons.numTriangles(); i < I; ++i) {
			const openvdb::Vec3I& quad = polygons.triangle(i);
			mesh.triIndexes.push_back(uint3(quad[2], quad[1], quad[0]));
		}
	}
	mesh.updateVertexNormals();
	mesh.updateBoundingBox();
}
bool WriteToRawFile(const std::string& fileName,openvdb::tools::Dense<float,openvdb::tools::MemoryLayout::LayoutZYX>& dense){
	std::ostringstream vstr;
	vstr << fileName << ".raw";
	FILE* f = fopen(vstr.str().c_str(), "wb");
	openvdb::CoordBBox bbox =dense.bbox();
	Coord dims = bbox.max() - bbox.min() + Coord(1, 1, 1);
	openvdb::Coord P(0, 0, 0);
	for (P[2] = bbox.min()[2]; P[2] <= bbox.max()[2]; ++P[2]) {
		for (P[1] = bbox.min()[1]; P[1] <= bbox.max()[1]; ++P[1]) {

			for (P[0] = bbox.min()[0]; P[0] <= bbox.max()[0]; ++P[0]) {
				float val = dense.getValue(P);
				fwrite(&val, sizeof(float), 1, f);
			}
		}
	}
	fclose(f);
	std::stringstream sstr;
	sstr << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	sstr << "<!-- MIPAV header file -->\n";
	sstr
			<< "<image xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" nDimensions=\"3\">\n";
	sstr << "	<Dataset-attributes>\n";
	sstr << "		<Image-offset>0</Image-offset>\n";
	sstr << "		<Data-type>Float</Data-type>\n";
	sstr << "		<Endianess>Little</Endianess>\n";
	sstr << "		<Extents>" << dims[0] << "</Extents>\n";
	sstr << "		<Extents>" << dims[1] << "</Extents>\n";
	sstr << "		<Extents>" << dims[2] << "</Extents>\n";
	sstr << "		<Resolutions>\n";
	sstr << "			<Resolution>1.0</Resolution>\n";
	sstr << "			<Resolution>1.0</Resolution>\n";
	sstr << "			<Resolution>1.0</Resolution>\n";
	sstr << "		</Resolutions>\n";
	sstr << "		<Slice-spacing>1.0</Slice-spacing>\n";
	sstr << "		<Slice-thickness>0.0</Slice-thickness>\n";
	sstr << "		<Units>Millimeters</Units>\n";
	sstr << "		<Units>Millimeters</Units>\n";
	sstr << "		<Units>Millimeters</Units>\n";
	sstr << "		<Compression>none</Compression>\n";
	sstr << "		<Orientation>Unknown</Orientation>\n";
	sstr << "		<Subject-axis-orientation>Unknown</Subject-axis-orientation>\n";
	sstr << "		<Subject-axis-orientation>Unknown</Subject-axis-orientation>\n";
	sstr << "		<Subject-axis-orientation>Unknown</Subject-axis-orientation>\n";
	sstr << "		<Origin>0.0</Origin>\n";
	sstr << "		<Origin>0.0</Origin>\n";
	sstr << "		<Origin>0.0</Origin>\n";
	sstr << "		<Modality>Unknown Modality</Modality>\n";
	sstr << "	</Dataset-attributes>\n";
	sstr << "</image>\n";
	std::ofstream myfile;
	std::stringstream xmlFile;
	xmlFile << fileName << ".xml";
	myfile.open(xmlFile.str().c_str(), std::ios_base::out);
	myfile << sstr.str();
	myfile.close();
	std::cout << xmlFile.str() << std::endl;
	return true;
}
bool WriteToRawFile(const std::string& fileName,openvdb::tools::Dense<openvdb::Vec3s,openvdb::tools::MemoryLayout::LayoutZYX>& dense){
	std::ostringstream vstr;
	vstr << fileName << ".raw";
	FILE* f = fopen(vstr.str().c_str(), "wb");
	openvdb::CoordBBox bbox =dense.bbox();
	Coord dims = bbox.max() - bbox.min() + Coord(1, 1, 1);
	openvdb::Coord P(0, 0, 0);
	for (P[2] = bbox.min()[2]; P[2] <= bbox.max()[2]; ++P[2]) {
		for (P[1] = bbox.min()[1]; P[1] <= bbox.max()[1]; ++P[1]) {
			for (P[0] = bbox.min()[0]; P[0] <= bbox.max()[0]; ++P[0]) {
				Vec3f val = dense.getValue(P);
				fwrite(&val[0], sizeof(float), 1, f);
			}
		}
	}

	for (P[2] = bbox.min()[2]; P[2] <= bbox.max()[2]; ++P[2]) {
		for (P[1] = bbox.min()[1]; P[1] <= bbox.max()[1]; ++P[1]) {
			for (P[0] = bbox.min()[0]; P[0] <= bbox.max()[0]; ++P[0]) {
				Vec3f val = dense.getValue(P);
				fwrite(&val[1], sizeof(float), 1, f);
			}
		}
	}

	for (P[2] = bbox.min()[2]; P[2] <= bbox.max()[2]; ++P[2]) {
		for (P[1] = bbox.min()[1]; P[1] <= bbox.max()[1]; ++P[1]) {
			for (P[0] = bbox.min()[0]; P[0] <= bbox.max()[0]; ++P[0]) {
				Vec3f val = dense.getValue(P);
				fwrite(&val[2], sizeof(float), 1, f);
			}
		}
	}

	fclose(f);
	std::stringstream sstr;
	sstr << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	sstr << "<!-- MIPAV header file -->\n";
	sstr
			<< "<image xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" nDimensions=\"4\">\n";
	sstr << "	<Dataset-attributes>\n";
	sstr << "		<Image-offset>0</Image-offset>\n";
	sstr << "		<Data-type>Float</Data-type>\n";
	sstr << "		<Endianess>Little</Endianess>\n";
	sstr << "		<Extents>" << dims[0] << "</Extents>\n";
	sstr << "		<Extents>" << dims[1] << "</Extents>\n";
	sstr << "		<Extents>" << dims[2] << "</Extents>\n";
	sstr << "		<Extents>3</Extents>\n";
	sstr << "		<Resolutions>\n";
	sstr << "			<Resolution>1.0</Resolution>\n";
	sstr << "			<Resolution>1.0</Resolution>\n";
	sstr << "			<Resolution>1.0</Resolution>\n";
	sstr << "		</Resolutions>\n";
	sstr << "		<Slice-spacing>1.0</Slice-spacing>\n";
	sstr << "		<Slice-thickness>0.0</Slice-thickness>\n";
	sstr << "		<Units>Millimeters</Units>\n";
	sstr << "		<Units>Millimeters</Units>\n";
	sstr << "		<Units>Millimeters</Units>\n";
	sstr << "		<Compression>none</Compression>\n";
	sstr << "		<Orientation>Unknown</Orientation>\n";
	sstr << "		<Subject-axis-orientation>Unknown</Subject-axis-orientation>\n";
	sstr << "		<Subject-axis-orientation>Unknown</Subject-axis-orientation>\n";
	sstr << "		<Subject-axis-orientation>Unknown</Subject-axis-orientation>\n";
	sstr << "		<Origin>0.0</Origin>\n";
	sstr << "		<Origin>0.0</Origin>\n";
	sstr << "		<Origin>0.0</Origin>\n";
	sstr << "		<Modality>Unknown Modality</Modality>\n";
	sstr << "	</Dataset-attributes>\n";
	sstr << "</image>\n";
	std::ofstream myfile;
	std::stringstream xmlFile;
	xmlFile << fileName << ".xml";
	myfile.open(xmlFile.str().c_str(), std::ios_base::out);
	myfile << sstr.str();
	myfile.close();
	std::cout << xmlFile.str() << std::endl;
	return true;
}
