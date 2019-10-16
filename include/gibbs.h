/*
 * gibbs.h
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

/// \file gibbs.h
/// \brief Gibbs sampling
/// \author Radu Marinescu radu.marinescu@ie.ibm.com
 
#ifndef IBM_MERLIN_GIBBS_H_
#define IBM_MERLIN_GIBBS_H_

#include "algorithm.h"
#include "graphical_model.h"

namespace merlin {

///
/// \brief Factor graph algorithm specialization for Gibbs sampling.
///
class gibbs: public graphical_model, public algorithm {
public:
	typedef graphical_model::findex findex;		///< Factor index
	typedef graphical_model::vindex vindex;		///< Variable index
	typedef graphical_model::flist flist; 		///< Collection of factor indices

public:

	///
	/// \brief Creates an empty Gibbs sampler.
	///
	gibbs() :
			graphical_model() {
		set_properties();
	}

	///
	/// \brief Creates a Gibbs sampler from an existing graphical model.
	///
	gibbs(const graphical_model& fg) :
			graphical_model(fg) {
		set_properties();
	}

	///
	/// \brief Clone the Gibbs sampler.
	/// \return the pointer to the cloned sampler.
	///
//	gibbs* clone() const {
//		gibbs* fg = new gibbs(*this);
//		return fg;
//	}

	///
	/// \brief Properties of the sampler.
	///
	MER_ENUM( Property , Task,TempMin,TempMax,Iter,Samples,Debug );

	///
	/// \brief Set the properties of the sampler.
	/// \param opt 	The string of comma separated property value pairs
	///
	virtual void set_properties(std::string opt = std::string()) {
		if (opt.length() == 0) {
			set_properties(
					"Task=MAR,TempMin=1.0,TempMax=1.0,Iter=10,Samples=100,Debug=0");
			m_order.clear();
			m_order.reserve(nvar());
			for (size_t v = 0; v < nvar(); v++)
				m_order.push_back(v);
			return;
		}

		m_debug = false;
		std::vector<std::string> strs = split(opt, ',');
		for (size_t i = 0; i < strs.size(); ++i) {
			std::vector<std::string> asgn = split(strs[i], '=');
			switch (Property(asgn[0].c_str())) {
			case Property::Task:
				m_task = Task(asgn[1].c_str());
				break;
			case Property::TempMin:
				m_temp_min = strtod(asgn[1].c_str(), NULL);
				break;
			case Property::TempMax:
				m_temp_max = strtod(asgn[1].c_str(), NULL);
				break;
			case Property::Iter:
				m_num_iter = atol(asgn[1].c_str());
				break;
			case Property::Samples:
				m_num_samples = atol(asgn[1].c_str());
				break;
			case Property::Debug:
				m_debug = (atol(asgn[1].c_str()) == 0 ? false : true);
				break;
			default:
				break;
			}
		}
	}

	///
	/// \brief Initialize the Gibbs sampler.
	///
	void init();

	///
	/// \brief Run the Gibbs sampler.
	///
	void run();

	double lb() const {
		return m_lb;
	}
	double ub() const {
		throw std::runtime_error("Not available");
	}
	std::vector<index> best_config() const {
		return m_best_config;
	}
	double logZ() const {
		throw std::runtime_error("Not available");
	}
	double logZub() const {
		throw std::runtime_error("Not available");
	}
	double logZlb() const {
		throw std::runtime_error("Not available");
	}

	///
	/// \brief Return the list of samples.
	///
	const std::vector<std::vector<index> >& samples() {
		return m_samples;
	}

	const factor& belief(size_t i) const {
		return m_beliefs[i];
	}
	const factor& belief(variable v) const {
		//return belief(with_variable(v)[0]);
		return m_beliefs[v];
	}
	const factor& belief(variable_set vs) const {
		throw std::runtime_error("Not implemented");
	}
	const std::vector<factor>& beliefs() const {
		return m_beliefs;
	}

	///
	/// \brief Write the solution to the output file.
	/// \param filename 	The output file name
	/// \param evidence 	The evidence variable value pairs
	/// \param old2new		The mapping between old and new variable indexing
	/// \param orig 		The graphical model prior to asserting evidence
	///
	void write_solution(const char* file_name, const std::map<size_t, size_t>& evidence,
			const std::map<size_t, size_t>& old2new, const graphical_model& orig );

	///
	/// \brief Write the solution to the output file.
	/// \param filename 	The output file name
	/// \param evidence 	The evidence variable value pairs
	/// \param old2new		The mapping between old and new variable indexing
	/// \param orig 		The graphical model prior to asserting evidence
	///
	void write_solution(std::ostream& out, const std::map<size_t, size_t>& evidence,
			const std::map<size_t, size_t>& old2new, const graphical_model& orig,
			const std::set<size_t>& dummies, int output_format);

	///
	/// \brief Inference tasks supported.
	///
	MER_ENUM( Task, PR,MAR,MAP );

private:
	// Members:
	Task m_task;							///< Inference task

	size_t m_num_samples;
	size_t m_num_iter;
	std::vector<index> m_state;
	variable_order_t m_order;
	std::vector<index> m_best_config;
	double m_lb;
	std::vector<factor> m_beliefs;
	std::vector<std::vector<index> > m_samples;
	double m_temp_min;
	double m_temp_max;
	double m_temp;
	bool m_debug;								///< Internal debugging flag
};

} // namespace

#endif // re-include

