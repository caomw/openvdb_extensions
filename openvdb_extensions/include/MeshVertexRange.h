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

#ifndef INCLUDE_MESHVERTEXRANGE_H_
#define INCLUDE_MESHVERTEXRANGE_H_
#include "Constellation.h"

class MeshVertexRange {
public:
	class Iterator {
	public:
		Iterator(const MeshVertexRange& range, size_t pos) :
				mRange(range), mPos(pos) {
			assert(this->isValid());
		}

		/// Advance to the next leaf node.
		Iterator& operator++() {
			++mPos;
			return *this;
		}
		/// Return a reference to the leaf node to which this iterator is pointing.
		size_t operator*() const {
			return mPos;
		}

		/// Return a pointer to the leaf node to which this iterator is pointing.
		size_t operator->() const {
			return *this;
		}

		/// Return the index into the leaf array of the current leaf node.
		size_t pos() const {
			return mPos;
		}
		bool isValid() const {
			return mPos >= mRange.mBegin && mPos <= mRange.mEnd;
		}
		/// Return @c true if this iterator is not yet exhausted.
		bool test() const {
			return mPos < mRange.mEnd;
		}
		/// Return @c true if this iterator is not yet exhausted.
		operator bool() const {
			return this->test();
		}
		/// Return @c true if this iterator is exhausted.
		bool empty() const {
			return !this->test();
		}
		bool operator!=(const Iterator& other) const {
			return (mPos != other.mPos) || (&mRange != &other.mRange);
		}
		bool operator==(const Iterator& other) const {
			return !(*this != other);
		}
		const MeshVertexRange& leafRange() const {
			return mRange;
		}

	protected:
		const MeshVertexRange& mRange;
		size_t mPos;
	}; // end Iterator

	MeshVertexRange(size_t begin, size_t end,Constellation& mesh,
			size_t grainSize = 1) :
				 mBegin(begin),mEnd(end), mGrainSize(grainSize), mMesh(
					mesh) {
	}

	MeshVertexRange(Constellation& mesh, size_t grainSize = 1) :
			 mBegin(0), mEnd(mesh.mVertexes.size()),mGrainSize(
					grainSize), mMesh(mesh) {
	}

	Iterator begin() const {
		return Iterator(*this, mBegin);
	}

	Iterator end() const {
		return Iterator(*this, mEnd);
	}

	size_t size() const {
		return mEnd - mBegin;
	}

	size_t grainsize() const {
		return mGrainSize;
	}

	const Constellation& getMesh() const {
		return mMesh;
	}

	bool empty() const {
		return !(mBegin < mEnd);
	}

	bool is_divisible() const {
		return mGrainSize < this->size();
	}

	MeshVertexRange(MeshVertexRange& r, tbb::split) :
		 mBegin(r.mBegin),mEnd(r.mEnd), mGrainSize(r.mGrainSize), mMesh(
					r.mMesh) {
		mEnd=r.mEnd;
		mBegin=doSplit(r);
	}
private:
	size_t mBegin, mEnd, mGrainSize;
	Constellation& mMesh;
	static size_t doSplit(MeshVertexRange& r) {
		assert(r.is_divisible());
		size_t middle = r.mBegin + (r.mEnd - r.mBegin) / 2u;
		r.mEnd = middle;
		return middle;
	}
};


#endif /* INCLUDE_MESHVERTEXRANGE_H_ */
