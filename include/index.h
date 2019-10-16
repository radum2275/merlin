/*
 * index.h
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

/// \file index.h
/// \brief Indexing routines
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_MERLIN_INDEX_H_
#define IBM_MERLIN_INDEX_H_

#include <iostream>

namespace merlin {

///
/// \brief Subindex for iterating over configurations of a set of variables.
///
class subindex {
public:
	typedef variable_set::vsize vsize;	///< Variable index

	vsize m_idx;				///< Current index position
	vsize m_end;				///< Last position
	vsize m_nd;					///< Number of variables in the sub-index
	vsize *m_state;				///< Vector of variable-indices (values) for the current position (1-based)
	const vsize *m_dims;		///< Dimensions of each variable
	bool *m_skipped;			///< Variables not included in the sub-index
	vsize *m_add;				///< How much to add when we increment
	vsize *m_subtract;			///< How much to subtract when we wrap each variable

	///
	/// \brief Construct the sub-index.
	///
	subindex(const variable_set& full, const variable_set& sub) {
		assert( full >> sub);
		m_idx = 0;
		m_end = 1;
		m_nd = full.nvar();
		m_dims = full.dims();
		m_state = new vsize[m_nd];
		m_add = new vsize[m_nd];
		m_subtract = new vsize[m_nd];
		m_skipped = new bool[m_nd];
		// Compute reference index updates
		vsize i, j;
		for (i = 0, j = 0; i < m_nd; ++i) {
			m_state[i] = 1;                       // start with [0] = (1,1,1,1...)
			m_skipped[i] = (j >= sub.nvar() || sub[j] != full[i]); // are we sub-indexing this variable?
			if (i == 0)
				m_add[i] = 1;              // how much does adding one to this var
			else
				m_add[i] = m_add[i - 1] * (m_skipped[i - 1] ? 1 : m_dims[i - 1]); //    add to our position?
			m_subtract[i] = m_add[i] * ((m_skipped[i] ? 1 : m_dims[i]) - 1); // how much does wrapping back to 1 remove from pos?
			if (!m_skipped[i])
				j++;
			m_end *= m_dims[i];
		}
		//end_--;
	}

	///
	/// \brief Copy (constructor) the sub-index.
	///
	subindex(const subindex& s) {
		m_dims = s.m_dims;
		m_idx = s.m_idx;
		m_end = s.m_end;
		m_nd = s.m_nd;
		m_state = new vsize[m_nd];
		std::copy(s.m_state, s.m_state + m_nd, m_state);
		m_add = new vsize[m_nd];
		std::copy(s.m_add, s.m_add + m_nd, m_add);
		m_subtract = new vsize[m_nd];
		std::copy(s.m_subtract, s.m_subtract + m_nd, m_subtract);
		m_skipped = new bool[m_nd];
		std::copy(s.m_skipped, s.m_skipped + m_nd, m_skipped);
	}

	///
	/// \brief Sub-index destructor.
	///
	~subindex(void) {
		delete[] m_state;
		delete[] m_skipped;
		delete[] m_add;
		delete[] m_subtract;
	}

	///
	/// \brief Reset the sub-index.
	///
	subindex& reset() {
		for (size_t i = 0; i < m_nd; ++i)
			m_state[i] = 1;
		m_idx = 0;
		return *this;
	}
	
	///
	/// \brief Return the last position in the sub-index.
	///
	size_t end() {
		return m_end;
	}

	///
	/// \brief Prefix addition operator.
	///
	subindex& operator++(void) {
		for (size_t i = 0; i < m_nd; ++i) { // for each variable
			if (m_state[i] == m_dims[i]) { // if we reached the maximum, wrap around to 1
				m_state[i] = 1;          // subtract wrap value from position
				if (!m_skipped[i])
					m_idx -= m_subtract[i]; // and continue to next variable in sequence
			} else {                    // otherwise, increment variable index
				++m_state[i];             // add to our current position
				if (!m_skipped[i])
					m_idx += m_add[i];   // and break (leave later vars the same)
				break;
			}
		}
		return *this;
	}

	/// 
	/// \brief Postfix addition using prefix and copy.
	///
	subindex operator++(int) {
		subindex s(*this);
		++(*this);
		return s;
	}

	///
	/// \brief Conversion to index value.
	///
	operator size_t() const {
		return m_idx;
	}

	///
	/// \brief Output operator.
	///
	friend std::ostream& operator<<(std::ostream &os, const subindex& s) {

		os << "Subindex:" << std::endl;
		os << s.m_idx << ", " << s.m_end << std::endl;
		os << s.m_nd << std::endl;
		for (size_t i = 0; i < s.m_nd; ++i) {
			os  << s.m_state[i] << " ";
		}
		os << std::endl;
		for (size_t i = 0; i < s.m_nd; ++i) {
			os << s.m_dims[i] << " ";
		}
		os << std::endl;
		for (size_t i = 0; i < s.m_nd; ++i) {
			os << s.m_skipped[i] << " ";
		}
		os << std::endl;
		for (size_t i = 0; i < s.m_nd; ++i) {
			os << s.m_add[i] << " ";
		}
		os << std::endl;
		for (size_t i = 0; i < s.m_nd; ++i) {
			os << s.m_subtract[i] << " ";
		}
		os << std::endl;
		return os;
	}

};

///
/// \brief Superindex for iterating over the configurations of a set of variables.
///
class superindex {
public:
	typedef variable_set::vsize vsize;	///< Variable index

	vsize m_idx;			///< Index
	vsize m_end;			///< End position
	vsize m_ns;			///< Number of states
	vsize *m_state;		///< vector of variable-indices (values) for the current position (1-based)
	const vsize *m_dims;	///< Dimensions of each variable
	vsize *m_add;		///< How much to add when we increment

	///
	/// \brief Construct the super index.
	///
	superindex(const variable_set& full, const variable_set& sub, size_t offset) {
		assert( full >> sub);
		m_idx = offset;
		m_end = full.num_states();
		m_ns = sub.nvar();
		const vsize* dimf = full.dims();
		m_dims = sub.dims();
		m_state = new vsize[m_ns];
		m_add = new vsize[m_ns];
		size_t i, j, d = 1;
		for (i = 0, j = 0; j < m_ns; ++i) {
			if (full[i] == sub[j]) {
				m_state[j] = 1;
				m_add[j] = d;
				j++;
			}
			d *= dimf[i];
		}
		m_end = m_add[m_ns - 1] * m_dims[m_ns - 1] + offset;
	}

	///
	/// \brief Copy (constructor) a super index.
	///
	superindex(const superindex& s) {
		m_dims = s.m_dims;
		m_idx = s.m_idx;
		m_end = s.m_end;
		m_ns = s.m_ns;
		m_state = new vsize[m_ns];
		std::copy(s.m_state, s.m_state + m_ns, m_state);
		m_add = new vsize[m_ns];
		std::copy(s.m_add, s.m_add + m_ns, m_add);
	}

	///
	/// \brief Super index destructor.
	///
	~superindex(void) {
		delete[] m_state;
		delete[] m_add;
	}

	///
	/// \brief Reset the super index.
	///
	superindex& reset() {
		for (size_t i = 0; i < m_ns; ++i)
			m_state[i] = 1;
		m_idx = 0;
		return *this;
	}

	///
	/// \brief Return the last position in the super index.
	///
	size_t end() {
		return m_end;
	}

	///
	/// \brief Prefix addition operator.
	///
	superindex& operator++(void) {
		for (size_t i = 0; i < m_ns; ++i) { // for each variable
			if (m_state[i] == m_dims[i] && i < m_ns - 1) { // if we reached the maximum, wrap around to 1
				m_state[i] = 1;             //
				m_idx -= m_add[i] * (m_dims[i] - 1);
			} else {                  // otherwise, increment variable index
				++m_state[i];           // add to our current position
				m_idx += m_add[i];
				break;
			}
		}
		return *this;
	}

	///
	/// \brief Postfix addition using prefix and copy.
	///
	superindex operator++(int) {
		superindex s(*this);
		++(*this);
		return s;
	}

	///
	/// \brief Conversion to index value.
	///
	operator size_t() const {
		return m_idx;
	}
};

///
/// \brief Permutation mapping from a variable set to an order.
///
class permute_index {
public:

	// Construct permutation mapping from VarSet -> Order
	// (bigEndian: is 1st variable largest stride?)
	///
	/// \brief Construct permutation mapping from VarSet -> Order.
	///
	permute_index(const std::vector<variable>& order, bool big_endian = false) {
		m_i = 0;
		variable_set _vs = variable_set(order.begin(), order.end(), order.size()); // compute implicit source order (VarSet)
		m_pi.resize(order.size());
		m_dim.resize(order.size());
		for (size_t j = 0; j < order.size(); ++j)
			m_dim[j] = _vs[j].states(); // save dimensions in source order (VarSet)
		for (size_t j = 0; j < order.size(); ++j) { // compute mapping from target order to source
			size_t jj = big_endian ? order.size() - 1 - j : j;
			for (size_t k = 0; k < order.size(); ++k)
				if (_vs[k] == order[j]) {
					m_pi[jj] = k;
					break;
				}
		}
	}

	///
	/// \brief Get target index corresponding to current or specified source index.
	///
	operator size_t() {
		return convert(m_i);
	}

	///
	/// \brief Set the index.
	///
	permute_index& set(size_t i) {
		m_i = i;
		return *this;
	}

	///
	/// \brief Convert a source index into a target index.
	///
	size_t convert(size_t i) {
		std::vector<size_t> I(m_dim.size());
		size_t r = 0, m = 1;
		for (size_t v = 0; v < m_dim.size(); ++v) {
			I[v] = i % m_dim[v];
			i -= I[v];
			i /= m_dim[v];
		}
		for (size_t j = 0; j < m_dim.size(); ++j) {
			r += m * I[m_pi[j]];
			m *= m_dim[m_pi[j]];
		}
		return r;
	}

	///
	/// \brief Invert mapping from Order -> VarSet.
	///
	permute_index inverse() {
		permute_index inv(*this);
		for (size_t i = 0; i < m_pi.size(); ++i) {
			inv.m_pi[m_pi[i]] = i;
			inv.m_dim[i] = m_dim[m_pi[i]];
		}
		inv.m_i = (size_t) *this;
		return inv;
	}

	// Can be used as an iterator as well:

	///
	/// \brief Iterate forward.
	///
	permute_index& operator++(void) {
		++m_i;
		return *this;
	}

	///
	/// \brief Iterate forward.
	///
	permute_index operator++(int) {
		permute_index r(*this);
		++m_i;
		return r;
	}

	///
	/// \brief Iterate backwards.
	///
	permute_index& operator--(void) {
		--m_i;
		return *this;
	}

	///
	/// \brief Iterate backwards.
	///
	permute_index operator--(int) {
		permute_index r(*this);
		--m_i;
		return r;
	}

private:
	size_t m_i;
	std::vector<size_t> m_pi;
	std::vector<size_t> m_dim;
};

///
/// \brief Convert index of variable configuration from order to set.
///
/// Given a configuration of the variables in the source order, convert the
/// corresponding index (0-based) to that corresponding to the configuration of
/// the same variables in a target order (i.e., permutation of the source vars).
///
class convert_index {
public:

	// BigEndian assumes that the first variable changes the fastest
	// UAI input is assumed to follow the LittleEndian convention, whereas
	// the internal representation of the factors assume BigEndian.

	///
	/// \brief Construct index convertor.
	/// \param order 	the source variable order
	/// \param source_big_endian the source configuration based on BigEndian
	/// \param target_big_endian the target configuration based on BigEndian
	///
	convert_index(const std::vector<variable>& order, bool source_big_endian = false,
			bool target_big_endian = false) {

		m_source_order = order; // save source variable order
		m_source_big_endian = source_big_endian; // source representation
		m_target_big_endian = target_big_endian; // target representation
		m_n = order.size(); // size of the source/target variable orders

		variable_set _vs = variable_set(order.begin(), order.end(), order.size());
		for (variable_set::const_iterator i = _vs.begin(); i != _vs.end(); ++i) {
			m_target_order.push_back(*i);
		}
	}

	convert_index(const std::vector<variable>& src_order, bool src_big_endian,
			const std::vector<variable>& trg_order, bool trg_big_endian) {
		m_source_order = src_order;
		m_target_order = trg_order;
		m_source_big_endian = src_big_endian;
		m_target_big_endian = trg_big_endian;
		m_n = src_order.size();
		assert(src_order.size() == trg_order.size());
	}

	///
	/// \brief Convert a source index into a target index.
	///
	size_t convert(size_t index) {

		std::map<size_t, size_t> config; // configuration of the source/target variables

		if (m_source_big_endian) { // start from the first variable in order
			for (size_t v = 0; v < m_n; ++v) {
				size_t k = m_source_order[v].states();
				size_t var = m_source_order[v].label();
				size_t val = index % k;
				index -= val;
				index /= k;
				config[var] = val;
			}
		} else { // start from last variable
			for (int v = m_n - 1; v >= 0; --v) {
				size_t k = m_source_order[v].states();
				size_t var = m_source_order[v].label();
				size_t val = index % k;
				index -= val;
				index /= k;
				config[var] = val;
			}
		}

		// compute the index of the configuration in the target space
		size_t r = 0, m = 1;
		if (m_target_big_endian) {
			for (size_t j = 0; j < m_n; ++j) {
				size_t var = m_target_order[j].label();
				size_t k = m_target_order[j].states();
				size_t val = config.find(var)->second;
				r += m*val;
				m *= k;
			}
		} else {
			for (int j = m_n - 1; j >= 0; --j) {
				size_t var = m_target_order[j].label();
				size_t k = m_target_order[j].states();
				size_t val = config.find(var)->second;
				r += m*val;
				m *= k;
			}
		}

		return r;
	}

private:
	size_t m_n;								///< Size of the source/target orders
	bool m_source_big_endian;				///< Source representation (BE or LE)
	bool m_target_big_endian;				///< Target representation (BE or LE)
	std::vector<variable> m_source_order;	///< Source variable order (unsorted)
	std::vector<variable> m_target_order;	///< Target variable set (sorted)
};

class index_config {

public:

	// Constructor
	index_config(const variable_set& vars, bool bigendian = true) {
		m_bigendian = bigendian;
		for (variable_set::const_iterator vi = vars.begin();
				vi != vars.end(); ++vi) {
			m_dims.push_back(vi->states());
			m_vars.push_back(*vi);
		}
	}

	// Get the configuration corresponding to an index
	std::map<size_t, size_t> convert(size_t index) {

		// Big endian: start from the first variable in order
		std::vector<size_t> I(m_dims.size()); // configuration of source variables
		if (m_bigendian) {
			for (size_t v = 0; v < m_dims.size(); ++v) {
				I[v] = index % m_dims[v];
				index -= I[v];
				index /= m_dims[v];
			}
		} else { // Little endian: start from the last variable in order
			for (int v = m_dims.size()-1; v >= 0; --v) {
				I[v] = index % m_dims[v];
				index -= I[v];
				index /= m_dims[v];
			}
		}

		std::map<size_t, size_t> config;
		for (size_t i = 0; i < I.size(); ++i) {
			size_t var = m_vars[i].label();
			size_t val = I[i];
			config.insert(std::make_pair(var, val));
		}

		return config;
	}

private:
	bool m_bigendian;						///< Big endian or little endian
	std::vector<size_t> m_dims;				///< Variable dimensions
	std::vector<variable> m_vars;			///< Source variable order (sorted)

};

class config_index {
public:
	config_index(const variable_set& vars, bool bigendian = true) {
		m_bigendian = bigendian;
		for (variable_set::const_iterator vi = vars.begin();
				vi != vars.end(); ++vi) {
			m_dims.push_back(vi->states());
			m_vars.push_back(*vi);
		}
	}

	size_t convert(std::map<size_t, size_t>& config) {
		// The source variables are a subset of the config vector

		std::vector<size_t> values;
		values.resize(m_vars.size());
		for (size_t i = 0; i < m_vars.size(); ++i) {
			size_t var = m_vars[i].label();
			std::map<size_t, size_t>::iterator mi = config.find(var);
			if (mi != config.end()) {
				size_t val = mi->second;
				values[i] = val;
			}
		}

		// Safety checks
		assert(values.size() == m_vars.size());

		size_t index = 0, offset = 1;
		if (m_bigendian) { // start from the first variable
			for (size_t i = 0; i < m_dims.size(); ++i) {
				index += values[i]*offset;
				offset *= m_dims[i];
			}
		} else {
			for (int i = m_dims.size()-1; i >= 0; --i) {
				index += values[i]*offset;
				offset *= m_dims[i];
			}
		}

		return index;
	}

private:
	bool m_bigendian;						///< Big endian or little endian
	std::vector<size_t> m_dims;				///< Variable dimensions
	std::vector<variable> m_vars;			///< Source variable order (sorted)
};

} // namespace

#endif /* IBM_MERLIN_INDEX_H_ */
