/*
 * set.h
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

/// \file set.h
/// \brief Set data structure
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_MERLIN_SET_H_
#define IBM_MERLIN_SET_H_

#include <assert.h>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <stdlib.h>
#include <stdint.h>

#include "vector.h"

namespace merlin {

///
/// \brief Set container storing unique elements following a specific order.
///
template<class T>
class my_set: protected my_vector<T> {
public:
	// Type definitions
	typedef typename my_vector<T>::iterator iterator;
	typedef typename my_vector<T>::const_iterator const_iterator;
	typedef typename my_vector<T>::reverse_iterator reverse_iterator;
	typedef typename my_vector<T>::const_reverse_iterator const_reverse_iterator;
	typedef std::out_of_range out_of_range;

	// Basic checks
	///
	/// \brief Return container size.
	///
	size_t size() const {
		return my_vector<T>::size();
	}

	///
	/// \brief Return maximum size of the container.
	///
	size_t max_size() const {
		return my_vector<T>::max_size();
	}

	///
	/// \brief Return the capacity of the container.
	///
	size_t capacity() const {
		return my_vector<T>::capacity();
	}

	///
	/// \brief Test whether container is empty.
	///
	bool empty() const {
		return my_vector<T>::empty();
	}

	// Setting, changing sizes

	///
	/// \brief Reserve memory for *n* elements.
	///
	void reserve(size_t n) {
		my_vector<T>::reserve(n);
	}

	///
	/// \brief Clear content.
	///
	void clear() {
		my_vector<T>::clear();
	}

	// Iterators

	///
	/// \brief Return const iterator to beginning.
	///
	const_iterator begin() const {
		return my_vector<T>::begin();
	}

	///
	/// \brief Return const iterator to end.
	///
	const_iterator end() const {
		return my_vector<T>::end();
	}

	///
	/// \brief Return const reverser iterator to reverse beginning.
	///
	const_reverse_iterator rbegin() const {
		return my_vector<T>::rbegin();
	}

	///
	/// \brief Return const reverse iterator to reverse end.
	///
	const_reverse_iterator rend() const {
		return my_vector<T>::rend();
	}

	// Accessor functions:

	///
	/// \brief Direct access to n-th element.
	///
	const T& operator[](size_t n) const {
		return my_vector<T>::operator[](n);
	}

	///
	/// \brief Direct access to n-th element.
	///
	const T& at(size_t n) const {
		return my_vector<T>::at(n);
	}

	///
	/// \brief Direct access to the last element.
	///
	const T& back() const {
		return my_vector<T>::back();
	}

	///
	/// \brief Direct access to the first element.
	///
	const T& front() const {
		return my_vector<T>::front();
	}

	/// 
	/// \brief Find a given element in the container.
	///
	const_iterator find(const T& x) const {
		return std::find(begin(), end(), x);
	}

	// Constructors  (missing a few from std::vector):

	///
	/// \brief Default constructor.
	///
	explicit my_set(size_t capacity = 10) : my_vector<T>() {
		reserve(capacity);
	};

	///
	/// \brief Copy-constructor.
	///
	my_set(const my_set& v) :	my_vector<T>() {
		my_vector<T>::operator=((my_vector<T>&) v);
	} 

	///
	/// \brief Constructor with input iterators.
	///
	template<class inIter>
	my_set(inIter first, inIter last) :
			my_vector<T>() {
		while (first != last) {
			*this |= *first;
			++first;
		}
	}

	///
	/// \brief Assignment operator.
	///
	my_set<T>& operator=(const my_set<T>& s) {
		my_vector<T>::operator=((my_vector<T>&) s);
		return *this;
	}

	///
	/// \brief Set destructor.
	///
	~my_set() {
	}

	// "Safe" insertion & removal
	my_set<T> operator|(const my_set<T>&) const;    ///< Union
	my_set<T> operator+(const my_set<T>&) const;    ///< Union
	my_set<T> operator/(const my_set<T>&) const;    ///< Set-diff
	my_set<T> operator-(const my_set<T>&) const;    ///< Set-diff
	my_set<T> operator&(const my_set<T>&) const;    ///< Intersect
	my_set<T> operator^(const my_set<T>&) const;    ///< Symm set diff
	my_set<T> operator|(const T&) const;         ///< Union
	my_set<T> operator+(const T&) const;    ///< Union
	my_set<T> operator/(const T&) const;    ///< Set-diff
	my_set<T> operator-(const T&) const;    ///< Set-diff
	my_set<T> operator&(const T&) const;    ///< Intersect
	my_set<T> operator^(const T&) const;    ///< Symm set diff

	my_set<T>& operator|=(const my_set<T>&);    ///< Union
	my_set<T>& operator+=(const my_set<T>&);    ///< Union
	my_set<T>& operator/=(const my_set<T>&);    ///< Set-diff
	my_set<T>& operator-=(const my_set<T>&);    ///< Set-diff
	my_set<T>& operator&=(const my_set<T>&);    ///< Intersect
	my_set<T>& operator^=(const my_set<T>&);    ///< Symm set diff

	my_set<T>& operator|=(const T&);    ///< Union
	my_set<T>& operator+=(const T&);    ///< Union
	my_set<T>& operator/=(const T&);    ///< Set-diff
	my_set<T>& operator-=(const T&);    ///< Set-diff
	my_set<T>& operator&=(const T&);    ///< Intersect
	my_set<T>& operator^=(const T&);    ///< Symm set diff

	///
	/// \brief Remove a single element.
	///
	void remove(const T& t);

	///
	/// \brief Add a single element.
	///
	void add(const T& t);

	///
	/// \brief Swap content.
	///
	void swap(my_set<T> &v) {
		my_vector<T>::swap((my_vector<T>&) v);
	}

	// Tests for equality and lexicographical order:

	///
	/// \brief Equality operator over lexicographic order.
	///
	bool operator==(const my_set<T>& t) const {
		return (this->size() == t.size())
				&& std::equal(this->begin(), this->end(), t.begin());
	}

	///
	/// \brief Not-equal operator over lexicographic order.
	///
	bool operator!=(const my_set<T>& t) const {
		return !(*this == t);
	}

	///
	/// \brief Less-than operator over lexicographic order.
	///	
	bool operator<(const my_set<T>& t) const {
		return std::lexicographical_compare(this->begin(), this->end(),
				t.begin(), t.end());
	}

	///
	/// \brief Less-or-equal-than operator over lexicographic order.
	///	
	bool operator<=(const my_set<T>& t) const {
		return (*this == t || *this < t);
	}

	///
	/// \brief Greater-than operator over lexicographic order.
	///	
	bool operator>(const my_set<T>& t) const {
		return !(*this <= t);
	}

	///
	/// \brief Greater-or-equal-than operator over lexicographic order.
	///	
	bool operator>=(const my_set<T>& t) const {
		return !(*this > t);
	}

protected:

	///
	/// \brief Update container.
	///
	void update(void) {
		my_vector<T>::update();
	}

	///
	/// \brief Resize container.
	///
	void resize(size_t n, const T& t = T()) {
		my_vector<T>::resize(n, t);
	}

	///
	/// \brief Return iterator to beginning.
	///
	iterator _begin() {
		return my_vector<T>::begin();
	}

	///
	/// \brief Return iterator to end.
	///
	iterator _end() {
		return my_vector<T>::end();
	}

	///
	/// \brief Return reverse iterator to reverse beginning.
	///
	reverse_iterator _rbegin() {
		return my_vector<T>::rbegin();
	}

	///
	/// \brief Return reverse iterator to reverse end.
	///
	reverse_iterator _rend() {
		return my_vector<T>::rend();
	}
};

template<class T>
my_set<T> my_set<T>::operator|(const my_set<T>& b) const {
	my_set<T> d;
	d.resize(size() + b.size());                         // reserve enough space
	typename my_set<T>::iterator dend;
	dend = std::set_union(begin(), end(), b.begin(), b.end(), d._begin()); // use stl set function
	d.m_n = dend - d.begin();
	d.update();
	return d;
}
template<class T> my_set<T> my_set<T>::operator+(const my_set<T>& b) const {
	return *this | b;
}
;
template<class T>
my_set<T> my_set<T>::operator/(const my_set<T>& b) const {
	my_set<T> d;
	d.resize(size());                                    // reserve enough space
	typename my_set<T>::iterator dend;
	dend = std::set_difference(begin(), end(), b.begin(), b.end(), d._begin()); // use stl set function
	d.m_n = dend - d.begin();
	d.update();
	return d;
}
template<class T> my_set<T> my_set<T>::operator-(const my_set<T>& b) const {
	return *this / b;
}
;
template<class T>
my_set<T> my_set<T>::operator&(const my_set<T>& b) const {
	my_set<T> d;
	d.resize(size());                                    // reserve enough space
	typename my_set<T>::iterator dend;
	dend = std::set_intersection(begin(), end(), b.begin(), b.end(),
			d._begin());           // use stl set function
	d.m_n = dend - d.begin();
	d.update();
	return d;
}
template<class T>
my_set<T> my_set<T>::operator^(const my_set<T>& b) const {
	my_set<T> d;
	d.resize(size() + b.size());                         // reserve enough space
	typename my_set<T>::iterator dend;
	dend = std::set_symmetric_difference(begin(), end(), b.begin(), b.end(),
			d._begin());   // use stl set function
	d.m_n = dend - d.begin();
	d.update();
	return d;
}

template<class T> my_set<T>& my_set<T>::operator|=(const my_set<T>& b) {
	*this = *this | b;
	return *this;
}
template<class T> my_set<T>& my_set<T>::operator+=(const my_set<T>& b) {
	*this = *this + b;
	return *this;
}
template<class T> my_set<T>& my_set<T>::operator/=(const my_set<T>& b) {
	*this = *this / b;
	return *this;
}
template<class T> my_set<T>& my_set<T>::operator-=(const my_set<T>& b) {
	*this = *this - b;
	return *this;
}
template<class T> my_set<T>& my_set<T>::operator&=(const my_set<T>& b) {
	*this = *this & b;
	return *this;
}
template<class T> my_set<T>& my_set<T>::operator^=(const my_set<T>& b) {
	*this = *this ^ b;
	return *this;
}

template<class T>
my_set<T> my_set<T>::operator|(const T& b) const {
	my_set<T> d;
	size_t n = size();
	d.resize(n + 1);                                // reserve enough space
	typename my_set<T>::iterator dend;
	dend = std::set_union(begin(), end(), &b, (&b) + 1, d._begin()); // use stl set function
	d.m_n = dend - d.begin();
	d.update();
	return d;
}
template<class T> my_set<T> my_set<T>::operator+(const T& b) const {
	return *this | b;
}
;
template<class T>
my_set<T> my_set<T>::operator/(const T& b) const {
	my_set<T> d;
	size_t n = size();
	d.resize(n);                                    // reserve enough space
	typename my_set<T>::iterator dend;
	dend = std::set_difference(begin(), end(), &b, (&b) + 1, d._begin()); // use stl set function
	d.m_n = dend - d.begin();
	d.update();
	return d;
}
template<class T> my_set<T> my_set<T>::operator-(const T& b) const {
	return *this / b;
}
;
template<class T>
my_set<T> my_set<T>::operator&(const T& b) const {
	my_set<T> d;
	size_t n = size();
	d.resize(n);                                    // reserve enough space
	typename my_set<T>::iterator dend;
	dend = std::set_intersection(begin(), end(), &b, (&b) + 1, d._begin()); // use stl set function
	d.m_n = dend - d.begin();
	d.update();
	return d;
}
template<class T>
my_set<T> my_set<T>::operator^(const T& b) const {
	my_set<T> d;
	size_t n = size();
	d.resize(n + 1);                                // reserve enough space
	typename my_set<T>::iterator dend;
	dend = std::set_symmetric_difference(begin(), end(), &b, (&b) + 1,
			d._begin());   // use stl set function
	d.m_n = dend - d.begin();
	d.update();
	return d;
}

template<class T> my_set<T>& my_set<T>::operator|=(const T& b) {
	*this = *this | b;
	return *this;
}
template<class T> my_set<T>& my_set<T>::operator+=(const T& b) {
	*this = *this + b;
	return *this;
}
template<class T> my_set<T>& my_set<T>::operator/=(const T& b) {
	*this = *this / b;
	return *this;
}
template<class T> my_set<T>& my_set<T>::operator-=(const T& b) {
	*this = *this - b;
	return *this;
}
template<class T> my_set<T>& my_set<T>::operator&=(const T& b) {
	*this = *this & b;
	return *this;
}
template<class T> my_set<T>& my_set<T>::operator^=(const T& b) {
	*this = *this ^ b;
	return *this;
}

template<class T> void my_set<T>::add(const T& t) {
	*this |= t;
}
template<class T> void my_set<T>::remove(const T& t) {
	*this /= t;
}

} // namespace
#endif  // re-include
