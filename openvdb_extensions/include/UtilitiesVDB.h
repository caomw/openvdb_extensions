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
#ifndef INCLUDE_UTILITIESVDB_H_
#define INCLUDE_UTILITIESVDB_H_

#include <AlloyMesh.h>
#include <AlloyVector.h>
#include <AlloyVolume.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <vector>
#include "RegularGrid.h"
#include "MACGrid.h"
#include "Constellation.h"
void ConvertLevelSetToVolume(const openvdb::FloatGrid::Ptr& grid,aly::Volume1f& volume);
void ConvertVectorFieldToVolume(const openvdb::VectorGrid::Ptr& grid,aly::Volume3f& volume);
void ConvertLevelSetToMesh(const openvdb::FloatGrid::Ptr& grid,aly::Mesh& mesh);
void ConvertMeshToLevelSet(const aly::Mesh& mesh,openvdb::FloatGrid::Ptr& grid);
void ConvertLevelSetToConstellation(const openvdb::FloatGrid::Ptr& grid,Constellation& mesh);
void ConvertLevelSetToConstellation(const openvdb::FloatGrid::Ptr& grid,openvdb::tools::VolumeToMesh& mesher,Constellation& mesh);
void ConvertLevelSetToBBoxTree(const openvdb::FloatGrid::Ptr& grid,aly::Mesh& mesh,int minDepth=0,int maxDepth=3);
template <class T> void Convert(const aly::Vector<T,2>& in,std::vector<openvdb::math::Vec2<T>>& out){
	size_t N=in.size();
	if(N>0){
		out.resize(N);
		std::memcpy(out.data(),in.ptr(),in.typeSize()*N);
	} else {
		out.clear();
	}
}
template <class T> void Convert(const aly::Vector<T,3>& in,std::vector<openvdb::math::Vec3<T>>& out){
	size_t N=in.size();
	if(N>0){
		out.resize(N);
		std::memcpy(out.data(),in.ptr(),in.typeSize()*N);
	} else {
		out.clear();
	}
}
template <class T> void Convert(const aly::Vector<T,4>& in,std::vector<openvdb::math::Vec4<T>>& out){
	size_t N=in.size();
	if(N>0){
		out.resize(N);
		std::memcpy(out.data(),in.ptr(),in.typeSize()*N);
	} else {
		out.clear();
	}
}
template <class T> void Convert(const std::vector<openvdb::math::Vec2<T>>& in,aly::Vector<T,2>& out){
	size_t N=in.size();
	if(N>0){
		out.resize(N);
		std::memcpy(out.ptr(),in.data(),out.typeSize()*N);
	} else {
		out.clear();
	}
}
template <class T> void Convert(const std::vector<openvdb::math::Vec3<T>>& in,aly::Vector<T,3>& out){
	size_t N=in.size();
	if(N>0){
		out.resize(N);
		std::memcpy(out.ptr(),in.data(),out.typeSize()*N);
	} else {
		out.clear();
	}
}
template <class T> void Convert(const std::vector<openvdb::math::Vec4<T>>& in,aly::Vector<T,4>& out){
	size_t N=in.size();
	if(N>0){
		out.resize(N);
		std::memcpy(out.ptr(),in.data(),out.typeSize()*N);
	} else {
		out.clear();
	}
}
bool WriteToRawFile(const std::string& file,const RegularGrid<float>& dense);
bool WriteToRawFile(const std::string& file,const MACGrid<float>& mac);
bool WriteToRawFile(const std::string& file,const openvdb::VectorGrid::Ptr& grid);
bool WriteToRawFile(const std::string& fileName,const openvdb::FloatGrid::Ptr& grid);
bool WriteToRawFile(const std::string& fileName,const openvdb::Int32Grid::Ptr& grid);
bool WriteToRawFile(const std::string& fileName,openvdb::tools::Dense<openvdb::Vec3s,openvdb::tools::MemoryLayout::LayoutZYX>& dense);
bool WriteToRawFile(const std::string& fileName,openvdb::tools::Dense<float,openvdb::tools::MemoryLayout::LayoutZYX>& grid);
#endif /* INCLUDE_UTILITIESVDB_H_ */
