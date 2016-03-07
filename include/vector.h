/*
 * vector.h
 *
 *  Created on: 24 Mar 2015
 *      Author: radu
 */

/// \file vector.h
/// \brief Vector data structure
/// \author Radu Marinescu


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
