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
#ifndef UPWINDGRADIENT_H_
#define VDBTOOLS_INCLUDE_ADVECTIONFORCE_H_
#include <openvdb/openvdb.h>
#include <openvdb/tools/GridOperators.h>
#include <openvdb/math/FiniteDifference.h>
#include <openvdb/math/Stencils.h>
#include <openvdb/math/Maps.h>
#include <openvdb/math/Transform.h>
#include <AlloyMath.h>
#include <memory>
using namespace openvdb;
using namespace openvdb::tools;
using namespace openvdb::math;
template<DScheme DiffScheme>
struct ISAdvectionForce {
	//static const DScheme FD = BIAS_SCHEME<DiffScheme>::FD;
	//static const DScheme BD = BIAS_SCHEME<DiffScheme>::BD;

	// random access version
	template<typename Accessor> static Vec3<typename Accessor::ValueType> result(
			const Accessor& grid, const Coord& ijk) {
		typedef typename Accessor::ValueType ValueType;
		typedef Vec3<ValueType> Vec3Type;

		/*
		 ValueType DXPLUS = D1<FD>::inX(grid, ijk);
		 ValueType DYPLUS = D1<FD>::inY(grid, ijk);
		 ValueType DZPLUS = D1<FD>::inZ(grid, ijk);

		 ValueType DXMINUS = D1<BD>::inX(grid, ijk);
		 ValueType DYMINUS = D1<BD>::inY(grid, ijk);
		 ValueType DZMINUS = D1<BD>::inZ(grid, ijk);
		 Vec3Type vec=Vec3Type(
		 std::max(static_cast<ValueType>(0.0), DXPLUS) + std::min(static_cast<ValueType>(0.0), DXMINUS),
		 std::max(static_cast<ValueType>(0.0), DYPLUS) + std::min(static_cast<ValueType>(0.0), DYMINUS),
		 std::max(static_cast<ValueType>(0.0), DZPLUS) + std::min(static_cast<ValueType>(0.0), DZMINUS));
		 */
		Vec3Type vec = Vec3Type(D1<DiffScheme>::inX(grid, ijk),
				D1<DiffScheme>::inY(grid, ijk), D1<DiffScheme>::inZ(grid, ijk));
		//vec.normalize(1E-6f);
		double scale = -aly::clamp(
				grid.getValue(ijk) / (double) (openvdb::LEVEL_SET_HALF_WIDTH),
				-1.0, 1.0);
		vec = scale * vec;
		return vec;
	}

	// stencil access version
	template<typename StencilT> static Vec3<typename StencilT::ValueType> result(
			const StencilT& stencil) {
		typedef typename StencilT::ValueType ValueType;
		typedef Vec3<ValueType> Vec3Type;
		Vec3Type vec = Vec3Type(D1<DiffScheme>::inX(stencil),
				D1<DiffScheme>::inY(stencil), D1<DiffScheme>::inZ(stencil));
		//vec.normalize(1E-6f);
		double scale = -clamp(
				stencil.template getValue<0, 0, 0>()
						/ (double) (openvdb::LEVEL_SET_HALF_WIDTH), -1.0, 1.0);

		vec = scale * vec;
		return vec;
	}
};
template<typename MapType, DScheme DiffScheme>
struct AdvectionForce {
	// random access version
	template<typename Accessor>
	static typename openvdb::math::internal::ReturnValue<Accessor>::Vec3Type result(
			const MapType& map, const Accessor& grid, const Coord& ijk) {
		typedef typename openvdb::math::internal::ReturnValue<Accessor>::Vec3Type Vec3Type;

		Vec3d iGradient(ISAdvectionForce<DiffScheme>::result(grid, ijk));
		return Vec3Type(map.applyIJT(iGradient, ijk.asVec3d()));
	}

	// stencil access version
	template<typename StencilT>
	static typename openvdb::math::internal::ReturnValue<StencilT>::Vec3Type result(
			const MapType& map, const StencilT& stencil) {
		typedef typename openvdb::math::internal::ReturnValue<StencilT>::Vec3Type Vec3Type;

		Vec3d iGradient(ISAdvectionForce<DiffScheme>::result(stencil));
		return Vec3Type(
				map.applyIJT(iGradient, stencil.getCenterCoord().asVec3d()));
	}
};

// Partial template specialization of AdvectionForce
// translation, any order
template<DScheme DiffScheme>
struct AdvectionForce<TranslationMap, DiffScheme> {
	// random access version
	template<typename Accessor>
	static typename openvdb::math::internal::ReturnValue<Accessor>::Vec3Type result(
			const TranslationMap&, const Accessor& grid, const Coord& ijk) {
		return ISAdvectionForce<DiffScheme>::result(grid, ijk);
	}

	// stencil access version
	template<typename StencilT>
	static typename openvdb::math::internal::ReturnValue<StencilT>::Vec3Type result(
			const TranslationMap&, const StencilT& stencil) {
		return ISAdvectionForce<DiffScheme>::result(stencil);
	}
};

/// Full template specialization of AdvectionForce
/// uniform scale, 2nd order
template<DScheme DiffScheme>
struct AdvectionForce<UniformScaleMap, DiffScheme> {
	// random access version
	template<typename Accessor>
	static typename openvdb::math::internal::ReturnValue<Accessor>::Vec3Type result(
			const UniformScaleMap& map, const Accessor& grid,
			const Coord& ijk) {
		typedef typename openvdb::math::internal::ReturnValue<Accessor>::ValueType ValueType;
		typedef typename openvdb::math::internal::ReturnValue<Accessor>::Vec3Type Vec3Type;

		Vec3Type iGradient(ISAdvectionForce<DiffScheme>::result(grid, ijk));
		ValueType inv2dx = ValueType(map.getInvTwiceScale()[0]);
		return iGradient * inv2dx;
	}

	// stencil access version
	template<typename StencilT>
	static typename openvdb::math::internal::ReturnValue<StencilT>::Vec3Type result(
			const UniformScaleMap& map, const StencilT& stencil) {
		typedef typename openvdb::math::internal::ReturnValue<StencilT>::ValueType ValueType;
		typedef typename openvdb::math::internal::ReturnValue<StencilT>::Vec3Type Vec3Type;

		Vec3Type iGradient(ISAdvectionForce<DiffScheme>::result(stencil));
		ValueType inv2dx = ValueType(map.getInvTwiceScale()[0]);
		return iGradient * inv2dx;
	}
};

/// Full template specialization of AdvectionForce
/// uniform scale translate, 2nd order
template<DScheme DiffScheme>
struct AdvectionForce<UniformScaleTranslateMap, DiffScheme> {
	// random access version
	template<typename Accessor>
	static typename openvdb::math::internal::ReturnValue<Accessor>::Vec3Type result(
			const UniformScaleTranslateMap& map, const Accessor& grid,
			const Coord& ijk) {
		typedef typename openvdb::math::internal::ReturnValue<Accessor>::ValueType ValueType;
		typedef typename openvdb::math::internal::ReturnValue<Accessor>::Vec3Type Vec3Type;

		Vec3Type iGradient(ISAdvectionForce<DiffScheme>::result(grid, ijk));
		ValueType inv2dx = ValueType(map.getInvTwiceScale()[0]);
		return iGradient * inv2dx;
	}

	// stencil access version
	template<typename StencilT>
	static typename openvdb::math::internal::ReturnValue<StencilT>::Vec3Type result(
			const UniformScaleTranslateMap& map, const StencilT& stencil) {
		typedef typename openvdb::math::internal::ReturnValue<StencilT>::ValueType ValueType;
		typedef typename openvdb::math::internal::ReturnValue<StencilT>::Vec3Type Vec3Type;

		Vec3Type iGradient(ISAdvectionForce<DiffScheme>::result(stencil));
		ValueType inv2dx = ValueType(map.getInvTwiceScale()[0]);
		return iGradient * inv2dx;
	}
};

/// Full template specialization of AdvectionForce
/// scale, 2nd order
template<DScheme DiffScheme>
struct AdvectionForce<ScaleMap, DiffScheme> {
	// random access version
	template<typename Accessor>
	static typename openvdb::math::internal::ReturnValue<Accessor>::Vec3Type result(
			const ScaleMap& map, const Accessor& grid, const Coord& ijk) {
		typedef typename openvdb::math::internal::ReturnValue<Accessor>::ValueType ValueType;
		typedef typename openvdb::math::internal::ReturnValue<Accessor>::Vec3Type Vec3Type;

		Vec3Type iGradient(ISAdvectionForce<DiffScheme>::result(grid, ijk));
		return Vec3Type(ValueType(iGradient[0] * map.getInvTwiceScale()[0]),
				ValueType(iGradient[1] * map.getInvTwiceScale()[1]),
				ValueType(iGradient[2] * map.getInvTwiceScale()[2]));
	}

	// stencil access version
	template<typename StencilT>
	static typename openvdb::math::internal::ReturnValue<StencilT>::Vec3Type result(
			const ScaleMap& map, const StencilT& stencil) {
		typedef typename openvdb::math::internal::ReturnValue<StencilT>::ValueType ValueType;
		typedef typename openvdb::math::internal::ReturnValue<StencilT>::Vec3Type Vec3Type;

		Vec3Type iGradient(ISAdvectionForce<DiffScheme>::result(stencil));
		return Vec3Type(ValueType(iGradient[0] * map.getInvTwiceScale()[0]),
				ValueType(iGradient[1] * map.getInvTwiceScale()[1]),
				ValueType(iGradient[2] * map.getInvTwiceScale()[2]));
	}
};

/// Full template specialization of AdvectionForce
/// scale translate, 2nd order
template<DScheme DiffScheme>
struct AdvectionForce<ScaleTranslateMap, DiffScheme> {
	// random access version
	template<typename Accessor>
	static typename openvdb::math::internal::ReturnValue<Accessor>::Vec3Type result(
			const ScaleTranslateMap& map, const Accessor& grid,
			const Coord& ijk) {
		typedef typename openvdb::math::internal::ReturnValue<Accessor>::ValueType ValueType;
		typedef typename openvdb::math::internal::ReturnValue<Accessor>::Vec3Type Vec3Type;

		Vec3Type iGradient(ISAdvectionForce<DiffScheme>::result(grid, ijk));
		return Vec3Type(ValueType(iGradient[0] * map.getInvTwiceScale()[0]),
				ValueType(iGradient[1] * map.getInvTwiceScale()[1]),
				ValueType(iGradient[2] * map.getInvTwiceScale()[2]));
	}

	// Stencil access version
	template<typename StencilT>
	static typename openvdb::math::internal::ReturnValue<StencilT>::Vec3Type result(
			const ScaleTranslateMap& map, const StencilT& stencil) {
		typedef typename openvdb::math::internal::ReturnValue<StencilT>::ValueType ValueType;
		typedef typename openvdb::math::internal::ReturnValue<StencilT>::Vec3Type Vec3Type;

		Vec3Type iGradient(ISAdvectionForce<DiffScheme>::result(stencil));
		return Vec3Type(ValueType(iGradient[0] * map.getInvTwiceScale()[0]),
				ValueType(iGradient[1] * map.getInvTwiceScale()[1]),
				ValueType(iGradient[2] * map.getInvTwiceScale()[2]));
	}
};
/// @brief Compute the gradient of a scalar grid.
template<typename InGridT, typename MaskGridType = typename openvdb::tools::gridop::ToMaskGrid<InGridT>::Type, typename InterruptT = util::NullInterrupter>
class AdvectionForceGrid {
public:
	typedef InGridT InGridType;
	typedef typename ScalarToVectorConverter<InGridT>::Type OutGridType;

	AdvectionForceGrid(const InGridT& grid, InterruptT* interrupt = NULL) :
			mInputGrid(grid), mInterrupt(interrupt), mMask(NULL) {
	}

	AdvectionForceGrid(const InGridT& grid, const MaskGridType& mask,
			InterruptT* interrupt = NULL) :
			mInputGrid(grid), mInterrupt(interrupt), mMask(&mask) {
	}

	typename OutGridType::Ptr process(bool threaded = true) {
		Functor functor(mInputGrid, mMask, threaded, mInterrupt);
		processTypedMap(mInputGrid.transform(), functor);
		if (functor.mOutputGrid)
			functor.mOutputGrid->setVectorType(VEC_COVARIANT);
		return functor.mOutputGrid;
	}

protected:
	struct Functor {
		Functor(const InGridT& grid, const MaskGridType* mask, bool threaded,
				InterruptT* interrupt) :
				mThreaded(threaded), mInputGrid(grid), mInterrupt(interrupt), mMask(
						mask) {
		}

		template<typename MapT>
		void operator()(const MapT& map) {
			typedef AdvectionForce<MapT, DScheme::CD_2ND> OpT;
			gridop::GridOperator<InGridType, MaskGridType, OutGridType, MapT,
					OpT, InterruptT> op(mInputGrid, mMask, map, mInterrupt);
			mOutputGrid = op.process(mThreaded); // cache the result
		}

		const bool mThreaded;
		const InGridT& mInputGrid;
		typename OutGridType::Ptr mOutputGrid;
		InterruptT* mInterrupt;
		const MaskGridType* mMask;
	}; // Private Functor

	const InGridT& mInputGrid;
	InterruptT* mInterrupt;
	const MaskGridType* mMask;
};
// end of AdvectionForce class

template<typename GridType, typename InterruptT> inline typename ScalarToVectorConverter<
		GridType>::Type::Ptr
advectionForce(const GridType& grid, bool threaded, InterruptT* interrupt);

template<typename GridType, typename MaskT, typename InterruptT> inline typename ScalarToVectorConverter<
		GridType>::Type::Ptr
advectionForce(const GridType& grid, const MaskT& mask, bool threaded,
		InterruptT* interrupt);

template<typename GridType> inline typename ScalarToVectorConverter<GridType>::Type::Ptr advectionForce(
		const GridType& grid, bool threaded = true) {
	return advectionForce<GridType, util::NullInterrupter>(grid, threaded, NULL);
}

template<typename GridType, typename MaskT> inline typename ScalarToVectorConverter<
		GridType>::Type::Ptr advectionForce(const GridType& grid,
		const MaskT& mask, bool threaded = true) {
	return advectionForce<GridType, MaskT, util::NullInterrupter>(grid, mask,
			threaded, NULL);
}
template<typename GridType, typename InterruptT> inline typename ScalarToVectorConverter<
		GridType>::Type::Ptr advectionForce(const GridType& grid, bool threaded,
		InterruptT* interrupt) {
	AdvectionForceGrid<GridType, typename gridop::ToMaskGrid<GridType>::Type,InterruptT> op(grid, interrupt);
	return op.process(threaded);
}

template<typename GridType, typename MaskT, typename InterruptT> inline typename ScalarToVectorConverter<
		GridType>::Type::Ptr advectionForce(const GridType& grid,
		const MaskT& mask, bool threaded, InterruptT* interrupt) {
	AdvectionForceGrid<GridType, MaskT, InterruptT> op(grid, mask, interrupt);
	return op.process(threaded);
}

#endif /* UPWINDGRADIENT_H_ */
