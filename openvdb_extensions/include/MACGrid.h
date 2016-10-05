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
#ifndef INCLUDE_MACGRID_H_
#define INCLUDE_MACGRID_H_
#include "RegularGrid.h"
template<typename ValueT> struct MACGrid {
protected:
	RegularGrid<ValueT> mX, mY, mZ;
	const size_t mRows;
	const size_t mCols;
	const size_t mSlices;
	const float mVoxelSize;
	openvdb::math::Transform::Ptr mTransform;
public:

	MACGrid(const openvdb::Coord& dim, float voxelSize,ValueT value=0.0) :
			mX(dim[0] + 1, dim[1], dim[2], voxelSize,value),
			mY(dim[0], dim[1] + 1, dim[2], voxelSize ,value),
			mZ(dim[0], dim[1], dim[2] + 1, voxelSize ,value),mRows(dim[0]),mCols(dim[1]),mSlices(dim[2]),mVoxelSize(voxelSize) {
		mTransform=openvdb::math::Transform::createLinearTransform(mVoxelSize);
	}
	MACGrid(int rows, int cols, int slices, float voxelSize,ValueT value=0.0) :
			mX(rows + 1, cols, slices, voxelSize,value),
			mY(rows, cols + 1, slices, voxelSize ,value),
			mZ(rows, cols, slices + 1, voxelSize ,value),mRows(rows),mCols(cols),mSlices(slices),mVoxelSize(voxelSize) {
		mTransform=openvdb::math::Transform::createLinearTransform(mVoxelSize);
	}
	RegularGrid<ValueT>& operator[](size_t i) {
		return (&mX)[i];
	}
	const RegularGrid<ValueT>& operator[](size_t i) const {
		return (&mX)[i];
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
	inline openvdb::math::Vec3<ValueT> interpolate(const openvdb::Vec3s& p) const {
		openvdb::math::Vec3<ValueT> u;
		float scale=1.0/mVoxelSize;
		u[0] = mX.interpolate(scale*p[0], scale*p[1]-0.5, scale*p[2]-0.5);
		u[1] = mY.interpolate(scale*p[0]-0.5, scale*p[1], scale*p[2]-0.5);
		u[2] = mZ.interpolate(scale*p[0]-0.5, scale*p[1]-0.5, scale*p[2]);
		return u;
	}
	inline openvdb::math::Vec3<ValueT> maxInterpolate(const openvdb::Vec3f& position,float radius){
		openvdb::math::Vec3<ValueT> values[15];
		const float n3=1.0f/std::sqrt(3.0f);
		MACGrid<ValueT>& grid=*this;
		values[0]=grid.interpolate(position);
		values[1]=grid.interpolate(position+openvdb::Vec3f(+radius,0,0));
		values[2]=grid.interpolate(position+openvdb::Vec3f(-radius,0,0));
		values[3]=grid.interpolate(position+openvdb::Vec3f(0,0,+radius));
		values[4]=grid.interpolate(position+openvdb::Vec3f(0,0,-radius));
		values[5]=grid.interpolate(position+openvdb::Vec3f(0,+radius,0));
		values[6]=grid.interpolate(position+openvdb::Vec3f(0,-radius,0));

		values[7 ]=grid.interpolate(position+openvdb::Vec3f(+radius*n3,+radius*n3,+radius*n3));
		values[8 ]=grid.interpolate(position+openvdb::Vec3f(+radius*n3,-radius*n3,+radius*n3));
		values[9 ]=grid.interpolate(position+openvdb::Vec3f(-radius*n3,+radius*n3,+radius*n3));
		values[10]=grid.interpolate(position+openvdb::Vec3f(-radius*n3,-radius*n3,+radius*n3));
		values[11]=grid.interpolate(position+openvdb::Vec3f(+radius*n3,+radius*n3,-radius*n3));
		values[12]=grid.interpolate(position+openvdb::Vec3f(+radius*n3,-radius*n3,-radius*n3));
		values[13]=grid.interpolate(position+openvdb::Vec3f(-radius*n3,+radius*n3,-radius*n3));
		values[14]=grid.interpolate(position+openvdb::Vec3f(-radius*n3,-radius*n3,-radius*n3));

		ValueT maxVal=std::numeric_limits<ValueT>::min();
		openvdb::math::Vec3<ValueT> final=values[0];
		for(int i=0;i<15;i++){
			ValueT lsqr=values[i].lengthSqr();
			if(lsqr>maxVal){
				maxVal=lsqr;
				final=values[i];
			}
		}
		return final;
	}
};
inline void SANITY_CHECK_MAC_GRID(){
	MACGrid<float> grid(128,128,128,1.0f);
}



#endif /* INCLUDE_MACGRID_H_ */
