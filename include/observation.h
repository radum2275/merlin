/*
 * observation.h
 *
 *  Created on: 16 Sep 2019
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

/// \file observation.h
/// \brief An observation or evidence variable in graphical models
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_MERLIN_OBSERVATION_H_
#define IBM_MERLIN_OBSERVATION_H_

#include "base.h"

namespace merlin {

class observation {

protected:
	int m_variable;
	int m_value;
	std::vector<double> m_likelihood;
	bool m_observed;
	bool m_virtual;

public:

	inline int var() {
		return m_variable;
	}

	inline int val() {
		return m_value;
	}

	inline bool is_observed() {
		return m_observed;
	}

	inline bool is_virtual() {
		return m_virtual;
	}

	inline std::vector<double> likelihood() {
		return m_likelihood;
	}

public:

	// Default constructor
	observation() : m_variable(MERLIN_UNKNOWN),
			m_value(MERLIN_UNKNOWN),
			m_observed(false),
			m_virtual(false) {};

	// Missing observation (value)
	observation(int var) : m_variable(var),
			m_value(-1),
			m_observed(false),
			m_virtual(false) {};

	// Regular observation (plain evidence)
	observation(int var, int val) : m_variable(var),
			m_value(val),
			m_observed(true),
			m_virtual(false) {};

	// Virtual observation (virtual/likelihood observation)
	observation(int var, std::vector<double>& likelihood) : m_variable(var),
			m_value(MERLIN_UNKNOWN),
			m_likelihood(likelihood),
			m_observed(false),
			m_virtual(true) {};

	~observation() {};
};

} // end namespace


#endif /* CORE_OBSERVATION_H_ */
