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
#ifndef INCLUDE_REGULARGRID_H_
#define INCLUDE_REGULARGRID_H_
#include <AlloyVector.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/Dense.h>
template<typename ValueT> class RegularGrid: public openvdb::tools::Dense<
		ValueT, openvdb::tools::MemoryLayout::LayoutZYX> {
private:
	ValueT* mPtr;
	size_t mStrideX;
	size_t mStrideY;
	size_t mRows;
	size_t mCols;
	size_t mSlices;
	float mVoxelSize;
	openvdb::BBoxd mBoundingBox;
	openvdb::math::Transform::Ptr mTransform;
public:
	RegularGrid(const openvdb::CoordBBox& boundingBox) :
		openvdb::tools::Dense<ValueT,openvdb::tools::MemoryLayout::LayoutZYX>(boundingBox) {

		openvdb::Coord minPt=boundingBox.min();
		openvdb::Coord maxPt=boundingBox.max();
		mBoundingBox=openvdb::BBoxd(openvdb::Vec3d(minPt[0],minPt[1],minPt[2]),openvdb::Vec3d(maxPt[0],maxPt[1],maxPt[2]));
		const openvdb::Coord dims=boundingBox.max()-boundingBox.min();
		mPtr = this->data();
		mStrideX = this->xStride();
		mStrideY = this->yStride();
		mRows = dims[0];
		mCols = dims[1];
		mSlices = dims[2];
		mVoxelSize=(boundingBox.max()-boundingBox.min())[0]/dims[0];
		mTransform=openvdb::math::Transform::createLinearTransform(mVoxelSize);

	}
	RegularGrid(const openvdb::Coord& dims, const openvdb::BBoxd& boundingBox,
			ValueT value=0.0) :openvdb::tools::Dense<ValueT,openvdb::tools::MemoryLayout::LayoutZYX>(dims, openvdb::Coord(0)),mBoundingBox(boundingBox) {
		this->fill(value);
		mPtr = this->data();
		mStrideX = this->xStride();
		mStrideY = this->yStride();
		mRows = dims[0];
		mCols = dims[1];
		mSlices = dims[2];
		//Assume isotropic voxels!
		mVoxelSize=(boundingBox.max()-boundingBox.min())[0]/dims[0];
		mTransform=openvdb::math::Transform::createLinearTransform(mVoxelSize);
	}
	RegularGrid(int rows, int cols, int slices,float voxelSize,ValueT value=0.0) :
			openvdb::tools::Dense<ValueT,
					openvdb::tools::MemoryLayout::LayoutZYX>(
					openvdb::Coord(rows, cols, slices), openvdb::Coord(0)),
					mStrideX(this->xStride()),
					mStrideY(this->yStride()),
					mRows(rows),
					mCols(cols),
					mSlices(slices),
					mVoxelSize(voxelSize),
					mBoundingBox(openvdb::Vec3d(0,0,0),openvdb::Vec3d(voxelSize*rows,voxelSize*cols,voxelSize*slices)) {
		mPtr = this->data();
		mTransform=openvdb::math::Transform::createLinearTransform(mVoxelSize);
		this->fill(value);
	}
	RegularGrid(const openvdb::Coord& dims,float voxelSize,ValueT value=0.0):RegularGrid(dims[0],dims[1],dims[2],voxelSize,value){

	}

	inline openvdb::math::Transform& transform() {
		return *mTransform;
	}
	inline openvdb::math::Transform::Ptr transformPtr() {
		return mTransform;
	}
	inline void setTrasnfrom(openvdb::math::Transform::Ptr transform){
		mTransform=transform;
	}
	ValueT& operator()(size_t i, size_t j, size_t k) {
		assert((i>=0&&i<mRows));
		assert((j>=0&&j<mCols));
		assert((k>=0&&k<mSlices));
		return mPtr[i * mStrideX + j * mStrideY + k];
	}
	const ValueT& operator()(size_t i, size_t j, size_t k) const {
		assert((i>=0&&i<mRows));
		assert((j>=0&&j<mCols));
		assert((k>=0&&k<mSlices));
		return mPtr[i * mStrideX + j * mStrideY + k];
	}
	ValueT& operator()(const openvdb::Coord& ijk) {
		assert((ijk[0]>=0&&ijk[0]<mRows));
		assert((ijk[1]>=0&&ijk[1]<mCols));
		assert((ijk[2]>=0&&ijk[2]<mSlices));
		return mPtr[ijk[0] * mStrideX + ijk[1] * mStrideY + ijk[2]];
	}
	const ValueT& operator()(const openvdb::Coord& ijk) const {
		assert((ijk[0]>=0&&ijk[0]<mRows));
		assert((ijk[1]>=0&&ijk[1]<mCols));
		assert((ijk[2]>=0&&ijk[2]<mSlices));
		return mPtr[ijk[0] * mStrideX + ijk[1] * mStrideY + ijk[2]];
	}
	const openvdb::BBoxd& getBoundingBox() const {
		return mBoundingBox;
	}
	inline const size_t size() const {
		return mRows * mCols * mSlices;
	}
	inline const size_t rows() const {
		return mRows;
	}
	inline const size_t cols() const {
		return mCols;
	}
	inline const size_t slices() const {
		return mSlices;
	}
	inline const float voxelSize() const {
		return mVoxelSize;
	}

	inline const openvdb::Coord dimensions() const {
		return openvdb::Coord(mRows,mCols,mSlices);
	}

	inline ValueT interpolate(float x, float y, float z) const {
		x = aly::clamp(x,0.0f,(float)mRows);
		y = aly::clamp(y,0.0f,(float)mCols);
		z = aly::clamp(z,0.0f,(float)mSlices);
		int i = std::min((int)x,(int)mRows-1);
		int j = std::min((int)y,(int)mCols-1);
		int k = std::min((int)z,(int)mSlices-1);
		const RegularGrid<ValueT>& q=*this;
		return	(k+1-z)*(((i+1-x)*q(i,j,k)+(x-i)*q(i+1,j,k))*(j+1-y) + ((i+1-x)*q(i,j+1,k)+(x-i)*q(i+1,j+1,k))*(y-j)) +
				(z-k)*(((i+1-x)*q(i,j,k+1)+(x-i)*q(i+1,j,k+1))*(j+1-y) + ((i+1-x)*q(i,j+1,k+1)+(x-i)*q(i+1,j+1,k+1))*(y-j));
	}
	inline ValueT interpolateWorld(float x, float y, float z) const {
		openvdb::Vec3d pt(x,y,z);
		if(!mBoundingBox.isInside(pt))return (openvdb::LEVEL_SET_HALF_WIDTH+1.0f);
		ValueT val=interpolate(pt-mBoundingBox.min());
		return val;
	}
	inline ValueT interpolate(const openvdb::Vec3d& pt) const{
		double x = aly::clamp(pt[0],0.0,(double)mRows);
		double y = aly::clamp(pt[1],0.0,(double)mCols);
		double z = aly::clamp(pt[2],0.0,(double)mSlices);
		int i = std::min((int)x,(int)mRows-1);
		int j = std::min((int)y,(int)mCols-1);
		int k = std::min((int)z,(int)mSlices-1);
		const RegularGrid<ValueT>& q=*this;
		return	(k+1-z)*(((i+1-x)*q(i,j,k)+(x-i)*q(i+1,j,k))*(j+1-y) + ((i+1-x)*q(i,j+1,k)+(x-i)*q(i+1,j+1,k))*(y-j)) +
				(z-k)*(((i+1-x)*q(i,j,k+1)+(x-i)*q(i+1,j,k+1))*(j+1-y) + ((i+1-x)*q(i,j+1,k+1)+(x-i)*q(i+1,j+1,k+1))*(y-j));
	}
	void copyTo(RegularGrid<ValueT>& out) {
		ValueT* src = this->data();
		ValueT* dest = out.data();
		memcpy(dest, src, sizeof(ValueT) * size());
	}
	void add(RegularGrid<ValueT>& out) {
		ValueT* src = this->data();
		ValueT* dest = out.data();
		size_t N = size();
		#pragma omp parallel for
		for (size_t n = 0; n < N; n++) {
			src[n] += dest[n];
		}
	}
	void set(const ValueT out) {
		size_t N = size();
		ValueT* src = this->data();
		#pragma omp parallel for
		for (size_t n = 0; n < N; n++) {
			src[n] = out;
		}
	}
	void subtract(RegularGrid<ValueT>& out) {
		ValueT* src = this->data();
		ValueT* dest = out.data();
		size_t N = size();
		#pragma omp parallel for
		for (size_t n = 0; n < N; n++) {
			src[n] -= dest[n];
		}
	}
	void subtractFrom(RegularGrid<ValueT>& out) {
		ValueT* src = this->data();
		ValueT* dest = out.data();
		size_t N = size();
		#pragma omp parallel for
		for (size_t n = 0; n < N; n++) {
			src[n] = dest[n] - src[n];
		}
	}
};

inline void SANITY_CHECK_REGULAR_GRID(){
	RegularGrid<float> grid(128,128,128,1.0f);
}




#endif /* INCLUDE_REGULARGRID_H_ */
