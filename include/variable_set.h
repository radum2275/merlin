/*
 * variable_set.h
 *
 *  Created on: Feb 8, 2013
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

/// \file variable_set.h
/// \brief A set of variables for graphical models
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_MERLIN_VARSET_H_
#define IBM_MERLIN_VARSET_H_

#include "variable.h"

namespace merlin {

typedef std::vector<size_t> variable_order_t;

///
/// Define a set of variables (ie, the scope of a factor).
///
class variable_set {
public:
	typedef size_t vindex;   		///< Variable IDs
	typedef size_t vsize;    		///< Dimension (cardinality) of variables

protected:
	// Members:

	std::vector<vindex> m_v;		///< Variable IDs (sorted)
	std::vector<vsize> m_dlocal;	///< Non-const version (equals d_) if we allocated the dimensions ourselves
	const vsize* m_d;				///< Dimensions of the variables

	///
	/// \brief Check for local storage of dimensions.
	///
	bool islocal() {
		return ( m_d == &m_dlocal[0] );
	}

	///
	/// \brief "Reserved memory" constructor.
	///
	explicit variable_set(const size_t nv) :
			m_v(nv), m_dlocal(nv) {
		m_d = &m_dlocal[0];
	}

public:

	///
	/// A const iterator for the variable set.
	///
	class const_iterator :
			public std::iterator<std::bidirectional_iterator_tag, variable, ptrdiff_t,
			const variable*, const variable&> {
	private:
		size_t m_i;
		variable_set const* m_vs;
		mutable variable m_tmp;
	public:
		const_iterator() : m_tmp() {
			m_vs = NULL;
			m_i = 0;
		}
		const_iterator(variable_set const* vs, size_t i) : m_tmp() {
			m_vs = vs;
			m_i = i;
		}
		const variable operator*() const {
			return variable(m_vs->m_v[m_i], m_vs->m_d[m_i]);
		}
		const variable* operator->() const {
			m_tmp=variable(m_vs->m_v[m_i],m_vs->m_d[m_i]);
			return &m_tmp;
		}  //!!! hacky
		const_iterator& operator++(void) {
			++m_i;
			return *this;
		}
		const_iterator operator++(int) {
			++(*this);
			return const_iterator(m_vs, m_i - 1);
		}
		const_iterator& operator--(void) {
			--m_i;
			return *this;
		}
		const_iterator operator--(int) {
			--(*this);
			return const_iterator(m_vs, m_i + 1);
		}
		bool operator==(const_iterator ci) {
			return (m_vs == ci.m_vs && m_i == ci.m_i);
		}
		bool operator!=(const_iterator ci) {
			return !(*this == ci);
		}
	};

	typedef const_iterator const_reverse_iterator; ///< A reverse const iterator.

	// Constructors:

	///
	/// \brief Default constructor.
	///
	variable_set(void) :
			m_v(), m_dlocal(), m_d(NULL) {
	}

	///
	/// \brief Copy constructor.
	///
	variable_set(const variable_set& vs) :
			m_v(vs.m_v), m_dlocal(vs.m_d, vs.m_d + vs.size()) {
		m_d = &m_dlocal[0];
	}

	///
	/// \brief Constructor with a single variable.
	///
	variable_set(const variable& v) :
			m_v(1, v.label()), m_dlocal(1, v.states()) {
		m_d = &m_dlocal[0];
	}

	///
	/// \brief Constructor with a pair of variables.
	///
	variable_set(const variable& v, const variable& v2) :
			m_v(2), m_dlocal(2) {
		if (v.label() < v2.label()) {
			m_v[0] = v.label();
			m_dlocal[0] = v.states();
			m_v[1] = v2.label();
			m_dlocal[1] = v2.states();
			m_d = &m_dlocal[0];
		} else if (v2.label() < v.label()) {
			m_v[1] = v.label();
			m_dlocal[1] = v.states();
			m_v[0] = v2.label();
			m_dlocal[0] = v2.states();
			m_d = &m_dlocal[0];
		} else {
			m_v[0] = v.label();
			m_dlocal[0] = v.states();
			m_v.resize(1);
			m_dlocal.resize(1);
			m_d = &m_dlocal[0];
		}
	}

	///
	/// \brief Constructor with input iterators.
	///
	template<typename Iter>
	variable_set(Iter B, Iter E, size_t size_hint = 0) :
			m_v(), m_dlocal() {
		m_v.reserve(size_hint);
		m_dlocal.reserve(size_hint);
		m_d = &m_dlocal[0];
		for (; B != E; ++B)
			*this |= *B;
	}

	///
	/// \brief Destructor.
	///
	~variable_set(void) {
	}

	///
	/// \brief Assignment operator.
	///
	variable_set& operator=(const variable_set& B) {
		m_v = B.m_v;
		m_dlocal = std::vector < vsize > (B.m_d, B.m_d + B.size());
		m_d = &m_dlocal[0];
		return *this;
	}

	///
	/// \brief Swap two variable sets.
	///
	void swap(variable_set& v) {
		m_v.swap(v.m_v);
		m_dlocal.swap(v.m_dlocal);
		std::swap(m_d, v.m_d);
	}

	void clear(void) {
		m_v.clear();
		m_dlocal.clear();
		m_d = NULL;
	}

	// Accessors:

	///
	/// \brief Return the size of the set.
	///
	size_t size() const {
		return m_v.size();
	}

	///
	/// \brief Return the number of variables.
	///
	size_t nvar() const {
		return m_v.size();
	}

	///
	/// \brief Return the dimensions of variables.
	///
	const vsize* dims() const {
		return m_d;
	}

	///
	/// \brief Return the cartesian product of the domains.
	///
	size_t num_states() const {
		size_t ns = 1;
		for (size_t i = 0; i < size(); ++i)
			ns *= m_d[i];
		return (ns == 0) ? 1 : ns;
	}

	///
	/// \brief Index operator.
	///
	const variable operator[](size_t c) const {
		return variable(m_v[c], m_d[c]);
	}

	// Tests for equality and lexicographical order:

	///
	/// \brief Compare sets of variables under lexicographical order (A == B).
	///
	bool operator==(const variable_set& B) const {
		return (size() == B.size())
				&& std::equal(m_v.begin(), m_v.end(), B.m_v.begin());
	}

	///
	/// \brief Compare sets of variables under lexicographical order (A != B).
	///	
	bool operator!=(const variable_set& B) const {
		return !(*this == B);
	}

	///
	/// \brief Compare sets of variables under lexicographical order (A < B).
	///
	bool operator<(const variable_set& t) const {
		return std::lexicographical_compare(this->begin(), this->end(),
				t.begin(), t.end());
	}

	///
	/// \brief Compare sets of variables under lexicographical order (A <= B).
	///	
	bool operator<=(const variable_set& t) const {
		return (*this == t || *this < t);
	}

	///
	/// \brief Compare sets of variables under lexicographical order (A > B).
	///	
	bool operator>(const variable_set& t) const {
		return !(*this <= t);
	}

	///
	/// \brief Compare sets of variables under lexicographical order (A >= B).
	///	
	bool operator>=(const variable_set& t) const {
		return !(*this > t);
	}

	///
	/// \brief Return the iterator to the beginning of the set.
	///
	const_iterator begin() const {
		return const_iterator(this, 0);
	}

	///
	/// \brief Return the iterator to the end of the set.
	///
	const_iterator end() const {
		return const_iterator(this, size());
	}

	///
	/// \brief Return the reverse iterator to the beginning of the set.
	///
	const_reverse_iterator rbegin() const {
		return const_iterator(this, size() - 1);
	}

	///
	/// \brief Return the reverse iterator to the end of the set.
	///
	const_reverse_iterator rend() const {
		return const_iterator(this, size_t(-1));
	}

	// Operators:  union (+,|), setdiff (-,/), intersection (&), xor (^)

	// Union (+, |)

	///
	/// \brief Union over sets (+=).
	/// \param B 	The input set of variables
	/// \return the reference to the modified *this* object representing
	///		the union between *this* and set B received as input.
	///
	variable_set& operator+=(const variable_set& B) {
		return (*this = *this + B);
	}

	///
	/// \brief Union over sets (|).
	/// \param B 	The input set of variables
	/// \return the reference to the modified *this* object representing
	///		the union between *this* and set B received as input.
	///
	variable_set operator|(const variable_set& B) const {
		return *this + B;
	}

	///
	/// \brief Union over sets (|=).
	/// \param B 	The input set of variables
	/// \return the reference to the modified *this* object representing
	///		the union between *this* and set B received as input.
	///	
	variable_set& operator|=(const variable_set& B) {
		return (*this = *this | B);
	}

	// Set-diff (-, /)
	///
	/// \brief Difference over sets (/).
	/// \param B 	The input set of variables
	/// \return the reference to the modified *this* object representing
	///		the set-difference between *this* and set B received as input.
	///	
	variable_set operator/(const variable_set& B) const {
		return *this - B;
	}
	///
	/// \brief Difference over sets (/=).
	/// \param B 	The input set of variables
	/// \return the reference to the modified *this* object representing
	///		the set-difference between *this* and set B received as input.
	///	
	variable_set& operator/=(const variable_set& B) {
		return (*this -= B);
	}

	// Intersection (&)
	///
	/// \brief Intersection over sets (&=).
	/// \param B 	The input set of variables
	/// \return the reference to the modified *this* object representing
	///		the intersection between *this* and set B received as input.
	///		
	variable_set& operator&=(const variable_set& B) {
		return (*this = (*this & B));
	}

    // Set-symmetric diff (^ xor)
	///
	/// \brief Symmetric Difference over sets (^ xor).
	/// \param B 	The input set of variables
	/// \return the reference to the modified *this* object representing
	///		the set-symmetric difference (xor) between *this* and 
	///		set B received as input.
	///	    
	variable_set& operator^=(const variable_set& B) {
		return (*this = *this ^ B);
	}

	// Operators on single-element arguments:

	///
	/// \brief Union with a single variable.
	///
	variable_set operator+(const variable& B) const {
		return *this + variable_set(B);
	}

	///
	/// \brief Union with a single variable.
	///
	variable_set& operator+=(const variable& B) {
		return *this += variable_set(B);
	}

	///
	/// \brief Union with a single variable.
	///	
	variable_set operator|(const variable& B) const {
		return *this | variable_set(B);
	}

	///
	/// \brief Union with a single variable.
	///	
	variable_set& operator|=(const variable& B) {
		return *this |= variable_set(B);
	}

	///
	/// \brief Set-difference with a single variable.
	///
	variable_set operator-(const variable& B) const {
		return *this - variable_set(B);
	}

	///
	/// \brief Set-difference with a single variable.	
	///
	variable_set& operator-=(const variable& B) {
		return *this -= variable_set(B);
	}

	///
	/// \brief Set-difference with a single variable.
	///	
	variable_set operator/(const variable& B) const {
		return *this / variable_set(B);
	}

	///
	/// \brief Set-difference with a single variable.
	///	
	variable_set& operator/=(const variable& B) {
		return *this /= variable_set(B);
	}

	///
	/// \brief Intersection with a single variable.
	///
	variable_set operator&(const variable& B) const {
		return *this & variable_set(B);
	}

	///
	/// \brief Intersection with a single variable.
	///
	variable_set& operator&=(const variable& B) {
		return *this &= variable_set(B);
	}
	
	///
	/// \brief Set-symmetric difference with a single variable.
	///
	variable_set operator^(const variable& B) const {
		return *this ^ variable_set(B);
	}
	
	///
	/// \brief Set-symmetric difference with a single variable.
	///
	variable_set& operator^=(const variable& B) {
		return *this ^= variable_set(B);
	}

	// Union
	///
	/// \brief Union over sets (+).
	/// \param B 	The input set of variables
	/// \return the set representing
	///		the union between *this* and set B received as input.
	///	
	variable_set operator+(const variable_set& B) const {
		variable_set dest(size() + B.size());
		size_t i, j, k;
		for (i = 0, j = 0, k = 0; i < size() && j < B.size();) {
			if (m_v[i] == B.m_v[j]) {
				dest.m_dlocal[k] = (m_d[i] == 0 ? B.m_d[j] : m_d[i]);
				dest.m_v[k++] = m_v[i++];
				j++;
			} else if (m_v[i] > B.m_v[j]) {
				dest.m_dlocal[k] = B.m_d[j];
				dest.m_v[k++] = B.m_v[j++];
			} else {
				dest.m_dlocal[k] = m_d[i];
				dest.m_v[k++] = m_v[i++];
			}
		}
		while (i < size()) {
			dest.m_dlocal[k] = m_d[i];
			dest.m_v[k++] = m_v[i++];
		}
		while (j < B.size()) {
			dest.m_dlocal[k] = B.m_d[j];
			dest.m_v[k++] = B.m_v[j++];
		}
		dest.m_v.resize(k);
		dest.m_dlocal.resize(k);
		return dest;
	}

	// Set-diff
	///
	/// \brief Difference over sets (-).
	/// \param B 	The input set of variables
	/// \return the set representing
	///		the set-difference between *this* and set B received as input.
	///		
	variable_set operator-(const variable_set& B) const {
		variable_set dest(size());
		size_t i, j, k;
		for (i = 0, j = 0, k = 0; i < size() && j < B.size();) {
			if (m_v[i] < B.m_v[j]) {
				dest.m_dlocal[k] = m_d[i];
				dest.m_v[k++] = m_v[i++];
			} else if (m_v[i] == B.m_v[j]) {
				i++;
				j++;
			} else {
				j++;
			}
		}
		while (i < size()) {
			dest.m_dlocal[k] = m_d[i];
			dest.m_v[k++] = m_v[i++];
		}
		dest.m_v.resize(k);
		dest.m_dlocal.resize(k);
		return dest;
	}

	///
	/// \brief Difference over sets (-).
	/// \param B 	The input set of variables
	/// \return the reference to the modified *this* object representing
	///		the set-difference between *this* and set B received as input.
	///		
	variable_set& operator-=(const variable_set& B) {
		variable_set& dest = *this;
		size_t i, j, k;
		for (i = 0, j = 0, k = 0; i < size() && j < B.size();) {
			if (m_v[i] < B.m_v[j]) {
				dest.m_dlocal[k] = m_d[i];
				dest.m_v[k++] = m_v[i++];
			} else if (m_v[i] == B.m_v[j]) {
				i++;
				j++;
			} else {
				j++;
			}
		}
		while (i < size()) {
			dest.m_dlocal[k] = m_d[i];
			dest.m_v[k++] = m_v[i++];
		}
		dest.m_v.resize(k);
		dest.m_dlocal.resize(k);
		return dest;
	}

	///
	/// \brief Intersection over sets (&).
	/// \param B 	The input set of variables
	/// \return the reference to the modified *this* object representing
	///		the intersection between *this* and set B received as input.
	///		
	variable_set operator&(const variable_set& B) const {
		variable_set dest(size() > B.size() ? B.size() : size());
		size_t i, j, k;
		for (i = 0, j = 0, k = 0; i < size() && j < B.size();) {
			if (m_v[i] < B.m_v[j]) {
				i++;
			} else if (m_v[i] == B.m_v[j]) {
				dest.m_dlocal[k] = (m_d[i] == 0 ? B.m_d[j] : m_d[i]);
				dest.m_v[k++] = m_v[i++];
				j++;
			} else {
				j++;
			}
		}
		dest.m_v.resize(k);
		dest.m_dlocal.resize(k);
		return dest;
	}

	// Set-symmetric diff (xor)
	///
	/// \brief Symmetric difference over sets (xor).
	/// \param B 	The input set of variables
	/// \return the set representing
	///		the set-symmetric difference between *this* and set B received as input.
	///			
	variable_set operator^(const variable_set& B) const {
		variable_set dest(size() + B.size());
		size_t i, j, k;
		for (i = 0, j = 0, k = 0; i < size() && j < B.size();) {
			if (m_v[i] == B.m_v[j]) {
				i++;
				j++;
			} else if (m_v[i] > B.m_v[j]) {
				dest.m_dlocal[k] = B.m_d[j];
				dest.m_v[k++] = B.m_v[j++];
			} else {
				dest.m_dlocal[k] = m_d[i];
				dest.m_v[k++] = m_v[i++];
			}
		}
		while (i < size()) {
			dest.m_dlocal[k] = m_d[i];
			dest.m_v[k++] = m_v[i++];
		}
		while (j < B.size()) {
			dest.m_dlocal[k] = B.m_d[j];
			dest.m_v[k++] = B.m_v[j++];
		}
		dest.m_v.resize(k);
		dest.m_dlocal.resize(k);
		return dest;
	}

	///
	/// \brief Add a new variable to the set.
	///
	variable_set& insert(const variable& v) {
		return *this += v;
	}
	///
	/// \brief Remove a variable from the set.
	///
	variable_set& erase(const variable& v) {
		return *this -= v;
	}

	///
	/// \brief Check if another set includes the current set *this*.
	///
	bool operator<<(const variable_set& S) const {
		return std::includes(S.m_v.begin(), S.m_v.end(), m_v.begin(), m_v.end());
	}

	///
	/// \brief Check if the current set *this* includes another set.
	///
	bool operator>>(const variable_set& S) const {
		return std::includes(m_v.begin(), m_v.end(), S.m_v.begin(), S.m_v.end());
	}

	///
	/// \brief Check if the current set *this* intersects another set.
	///
	bool intersects(const variable_set& S) const {
		return (*this & S).size() > 0;
	}

	///
	/// \brief Check if the current set contains a given variable.
	///
	bool contains(const variable& v) const {
		return std::binary_search(m_v.begin(), m_v.end(), v.label());
	}

	///
	/// \brief Output operator.
	///
	friend std::ostream& operator<<(std::ostream &os, variable_set const& v) {
		os << "[";
		if (v.nvar() != 0) {
			os << v.m_v[0];
			for (size_t i = 1; i < v.nvar(); i++)
				os << "," << v.m_v[i];
		}
		os << "]";
		return os;
	}

};

template <class MapType>
inline size_t sub2ind(const variable_set& vs, const MapType& val) {
  if (vs.num_states()==0) throw std::runtime_error("sub2ind with uninitialized dimensions");
  size_t i=0,m=1;
  for (size_t v=0;v<vs.size();++v) {
    i += m*(size_t)val[vs[v]]; m*=vs[v].states();
  }
  return i;
}
template <class MapType>
inline size_t sub2ind(const variable_set& vs, MapType& val) {
  if (vs.num_states()==0) throw std::runtime_error("sub2ind with uninitialized dimensions");
  size_t i=0,m=1;
  for (size_t v=0;v<vs.size();++v) {
    i += m*(size_t)val[vs[v]]; m*=vs[v].states();
  }
  return i;
}
template <class MapType>
inline void ind2sub(const variable_set& vs, size_t i, MapType& val) {
  if (vs.num_states()==0) throw std::runtime_error("ind2sub with uninitialized dimensions");
  for (size_t v=0; v<vs.size(); ++v) {
    size_t rem = i%vs[v].states();
    val[vs[v]] = rem; i-=rem; i/=vs[v].states();
  }
  return val;
}

inline std::vector<size_t> ind2sub(const variable_set& vs, size_t i) {
  if (vs.num_states()==0) throw std::runtime_error("ind2sub with uninitialized dimensions");
  std::vector<size_t> val; val.resize(vs.size());
  for (size_t v=0; v<vs.size(); ++v) {
    val[v] = i%vs[v].states(); i-=val[v]; i/=vs[v].states();
  }
  return val;
}


} // namespace

#endif // re-include
