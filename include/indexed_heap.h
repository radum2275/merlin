/*
 * indexed_heap.h
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

/// \file indexed_heap.h
/// \brief Indexed heap data structure
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_MERLIN_INDEXEDHEAP_H_
#define IBM_MERLIN_INDEXEDHEAP_H_

#include <stdexcept>
#include "assert.h"

#include "base.h"

namespace merlin {

/**
 * Indexed Heap : a "reversible" heap from key (double) to unique values (uint, 0..N-1; some can be missing) 
 *   for which the value can also be used to look up (access) or remove an entry
 *
 * Mainly this data structure is used in factor- or edge- priority scheduling; the factor or edge index
 *   is used as the unique value and its priority is the key.  This gives a lightweight priority queue
 *   that still enables edges to be easily "reprioritized".
 */

class indexed_heap {

private:
	std::vector<double> m_p;		///< Store the key (a double), typically a priority
	std::vector<size_t> m_id;	///< Store the identifying value (uint) of this key
	std::vector<size_t> m_rev;	///< Reverse lookup from value to position in the heap

public:

	///
	/// \brief Default constructor.
	///
	indexed_heap() :
			m_p(), m_id(), m_rev() {
	}
	;

	///
	/// \brief Constructor.
	///
	template<class InputPIterator, class InputIIterator>
	indexed_heap(InputPIterator p_begin, InputPIterator p_end,
			InputIIterator i_begin) :
			m_p(), m_id(), m_rev() {
		size_t dist = distance(p_begin, p_end);
		size_t mxi;
		m_p.resize(dist);
		m_id.resize(dist);
		for (size_t i = 0; i < m_p.size(); ++i, ++p_begin, ++i_begin) {
			m_p[i] = *p_begin;
			m_id[i] = *i_begin;
			mxi = std::max(mxi, *i_begin);
		}
		m_rev.resize(mxi + 1);
		size_t ii = 1;
		for (std::vector<size_t>::iterator i = m_id.begin(); i != m_id.end(); ++i)
			m_rev[*i] = ii;

		for (size_t i = size() / 2; i != 0; --i)
			max_heapify(i);
	}

	///
	/// \brief Clear the heap.
	///
	void clear() {
		m_p.clear();
		m_id.clear();
		m_rev.clear();
	}

	///
	/// \brief Insert an element into the heap.
	/// \param p 	The value of the key
	/// \param r 	The id associated with the key
	///
	void insert(double p, size_t r) {
		if (m_rev.size() <= r)
			m_rev.resize(r + 1, 0);
		size_t i;
		if (m_rev[r] != 0) {
			i = m_rev[r];
			m_p[i - 1] = p;
			max_heapify(i);
		} else {
			m_p.push_back(p);
			m_id.push_back(r);
			m_rev[r] = size();
			i = size();
		}
		for (;;) {
			size_t parent = i / 2;
			if (parent > 0 && m_p[parent - 1] < m_p[i - 1]) {
				heap_swap(parent, i);
				i = parent;
			} else
				return;
		}
	}

	///
	/// \brief Remove an id from the heap.
	/// \param r 	The id to be removed
	///
	void erase(size_t r) {
		size_t I = m_rev[r];
		if (I == 0)
			return;		// already gone
		if (I == size()) {
			m_rev[r] = 0;
			m_p.pop_back();
			m_id.pop_back();
		} else {
			heap_swap(I, size());
			m_rev[r] = 0;
			m_p.pop_back();
			m_id.pop_back();
			max_heapify(I);
		}
	}

	///
	/// \brief Internal debugging method.
	///
	void debug() {
		std::cout << "P: ";
		for (size_t i = 0; i < m_p.size(); ++i)
			std::cout << m_p[i] << " ";
		std::cout << "\n";
		std::cout << "I: ";
		for (size_t i = 0; i < m_id.size(); ++i)
			std::cout << m_id[i] << " ";
		std::cout << "\n";
		std::cout << "R: ";
		for (size_t i = 0; i < m_rev.size(); ++i)
			std::cout << m_rev[i] << " ";
		std::cout << "\n";
	}

	///
	/// \brief Remove the top element from the heap.
	///
	void pop() {
		assert(size() > 0);
		erase(m_id[0]);
	}

	///
	/// \brief Return the top element from the heap (highest priority).
	///
	std::pair<double, size_t> top() {
		assert(size() > 0);
		return std::pair<double, size_t>(m_p[0], m_id[0]);
	}

	///
	/// \brief Return the size of the heap.
	///
	size_t size() {
		return m_p.size();
	}

	///
	/// \brief Check if the heap is empty.
	///
	bool empty() {
		return m_p.empty();
	}

private:

	///
	/// \brief Swap two elements in the heap.
	///
	void heap_swap(size_t i, size_t j) {
		std::swap(m_p[i - 1], m_p[j - 1]);
		std::swap(m_id[i - 1], m_id[j - 1]);
		std::swap(m_rev[m_id[i - 1]], m_rev[m_id[j - 1]]);
	}

	///
	/// \brief Make the heap.
	///
	void max_heapify(size_t i) {
		for (;;) {
			size_t left = 2 * i, right = 2 * i + 1, largest = i;
			if (left <= m_p.size() && m_p[left - 1] > m_p[largest - 1])
				largest = left;
			if (right <= m_p.size() && m_p[right - 1] > m_p[largest - 1])
				largest = right;
			if (largest == i)
				return;
			heap_swap(largest, i);
			i = largest;
		}
	}

};

} // namespace

#endif
