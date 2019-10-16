/*
 * vector.h
 *
 *  Created on: 24 Mar 2015
 *      Author: radu
 *
 * Copyright (c) 2015, International Business Machines Corporation
 * and University of California Irvine. All rights reserved.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/// \file vector.h
/// \brief Vector data structure
/// \author Radu Marinescu radu.marinescu@ie.ibm.com


#ifndef IBM_MERLIN_VECTOR_H_
#define IBM_MERLIN_VECTOR_H_

#include <assert.h>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <stdlib.h>
#include <stdint.h>

#include <vector>


namespace merlin {

///
/// \brief Container representing an array that can change in size.
///
template<class T>
class vector: public std::vector<T> {
public:
	
	// Typedefs
	typedef typename std::vector<T>::iterator iterator;
	typedef typename std::vector<T>::const_iterator const_iterator;
	typedef typename std::vector<T>::reverse_iterator reverse_iterator;
	typedef typename std::vector<T>::const_reverse_iterator const_reverse_iterator;

	///
	/// \brief Construct empty vector.
	///
	explicit vector() :	std::vector<T>() {
	}

	///
	/// \brief Construct vector with a given nuber of elements.
	///
	explicit vector(size_t n, const T& t = T()) : std::vector<T>(n, t) {
	}

	///
	/// \brief Copy constructor.
	/// \param v 	A vector object of the same type.
	///
	vector(vector<T> const& v) : std::vector<T>((std::vector<T> const&) v) {
	}

	///
	/// \brief Construct vector from input iterators.
	///
	template<class inIter> vector(inIter first, inIter last) :
			std::vector<T>(first, last) {
	}

	///
	/// \brief Assign content.
	/// \param v 	A vector object of the same type.
	///
	vector<T>& operator=(const vector<T>& v) {
		std::vector<T>::operator=((std::vector<T>&) v);
		return *this;
	}

	// Tests for equality and lexicographical order

	///
	/// \brief Equality operator under lexicographical order.
	///
	bool operator==(const vector<T>& t) const {
		return (this->size() == t.size())
				&& std::equal(this->begin(), this->end(), t.begin());
	}

	///
	/// \brief Not-equal operator under lexicographical order.
	///	
	bool operator!=(const vector<T>& t) const {
		return !(*this == t);
	}

	///
	/// \brief Less-than operator under lexicographical order.
	///	
	bool operator<(const vector<T>& t) const {
		return std::lexicographical_compare(this->begin(), this->end(),
				t.begin(), t.end());
	}

	///
	/// \brief Less-or-equal-than operator under lexicographical order.
	///	
	bool operator<=(const vector<T>& t) const {
		return (*this == t || *this < t);
	}

	///
	/// \brief Greater-than operator under lexicographical order.
	///	
	bool operator>(const vector<T>& t) const {
		return !(*this <= t);
	}

	///
	/// \brief Greater-or-equal-than operator under lexicographical order.
	///	
	bool operator>=(const vector<T>& t) const {
		return !(*this > t);
	}

protected:
	size_t m_n;	///< Number of elements.

	///
	/// \brief Update the size of container.
	///
	void update(void) {
		if (m_n)
			this->resize(m_n);
		else
			this->clear();
	}
};

} // namespace

#endif  // re-include
