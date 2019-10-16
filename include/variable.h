/*
 * variable.h
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

/// \file variable.h
/// \brief A variable for graphical models
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_MERLIN_VAR_H_
#define IBM_MERLIN_VAR_H_

#include "base.h"

namespace merlin {

///
/// \brief Define a variable for graphical models.
///
class variable {
private:

	size_t m_label;		///< Label of the variable (its unique ID - 0, 1, 2, 3 ...)
	size_t m_states;		///< Number of possible values (its domain)

public:

	///
	/// \brief Default constructor.
	///
	/// Create a variable with label 0 and 0 states.
	///
	variable() :
			m_label(0), m_states(0) {
	}

	///
	/// \brief Constructor.
	///
	/// Create a variable with a given label and number of states.
	///
	variable(size_t label, size_t states) :
			m_label(label), m_states(states) {
	}

	///
	/// \brief Return the label of the variable.
	///
	inline size_t label() const {
		return m_label;
	}

	///
	/// \brief Return the domains size of the variable.
	///
	inline size_t states() const {
		return m_states;
	}

	// Comparison operators (assume lexicographical order):

	///
	/// \brief Compare variables by their labels under a lexicographic order (a < b).
	///
	inline bool operator<(const variable& v) const {
		return (m_label < v.m_label);
	}

	///
	/// \brief Compare variables by their labels under a lexicographic order (a <= b).
	///	
	inline bool operator<=(const variable& v) const {
		return (m_label <= v.m_label);
	}

	///
	/// \brief Compare variables by their labels under a lexicographic order (a > b).
	///
	inline bool operator>(const variable& v) const {
		return (m_label > v.m_label);
	}

	///
	/// \brief Compare variables by their labels under a lexicographic order (a >= b).
	///
	inline bool operator>=(const variable& v) const {
		return (m_label >= v.m_label);
	}

	///
	/// \brief Compare variables by their labels under a lexicographic order (a == b).
	///	
	inline bool operator==(const variable& v) const {
		return (m_label == v.m_label);
	}

	///
	/// \brief Compare variables by their labels under a lexicographic order (a != b).
	///
	inline bool operator!=(const variable& v) const {
		return (m_label != v.m_label);
	}

	///
	/// \brief Index operator.
	///
	inline operator size_t() const {
		return m_label;
	}
};

} // namespace

#endif // re-include
