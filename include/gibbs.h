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
/// \author Radu Marinescu
 
#ifndef IBM_MERLIN_GIBBS_H_
#define IBM_MERLIN_GIBBS_H_

#include "algorithm.h"
#include "graphical_model.h"

namespace merlin {

///
/// \brief FactorGraph algorithm specialization for Gibbs sampling.
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
	gibbs* clone() const {
		gibbs* fg = new gibbs(*this);
		return fg;
	}

	///
	/// \brief Properties of the sampler.
	///
	MER_ENUM( Property , TempMin,TempMax,Best,Beliefs,nIter,nSamples );

	///
	/// \brief Set the properties of the sampler.
	/// \param opt 	The string of comma separated property value pairs
	///
	virtual void set_properties(std::string opt = std::string()) {
		if (opt.length() == 0) {
			set_properties(
					"TempMin=1.0,TempMax=1.0,Best=0,Beliefs=0,nIter=1000,nSamples=100");
			m_order.clear();
			m_order.reserve(nvar());
			for (size_t v = 0; v < nvar(); v++)
				m_order.push_back(v);
			return;
		}

		std::vector<std::string> strs = split(opt, ',');
		for (size_t i = 0; i < strs.size(); ++i) {
			std::vector<std::string> asgn = split(strs[i], '=');
			switch (Property(asgn[0].c_str())) {
			case Property::TempMin:
				m_temp_min = strtod(asgn[1].c_str(), NULL);
				break;
			case Property::TempMax:
				m_temp_max = strtod(asgn[1].c_str(), NULL);
				break;
			case Property::Best:
				m_no_best = !atol(asgn[1].c_str());
				break;
			case Property::Beliefs:
				m_no_beliefs = !atol(asgn[1].c_str());
				break;
			case Property::nIter:
				m_num_iter = atol(asgn[1].c_str());
				break;
			case Property::nSamples:
				m_num_samples = atol(asgn[1].c_str());
				break;
			default:
				break;
			}
		}
	}

	///
	/// \brief Initialize the Gibbs sampler.
	///
	void init() {
		m_samples.clear();
		if (m_init_state.size() > 0)
			m_state = m_init_state;
		else {
			m_state.resize(nvar());// might want to search for a good initialization
			for (variable_order_t::iterator i = m_order.begin(); i != m_order.end(); ++i)
				m_state[var(*i)] = randi(var(*i).states());
		}
		if (!m_no_beliefs) {
			m_beliefs.resize(num_factors());
			for (size_t f = 0; f < num_factors(); ++f)
				m_beliefs[f] = factor(get_factor(f).vars(), 0.0);
		}
		if (!m_no_best) {
			m_best_config = m_state;
			m_lb = 0.0;
			for (size_t f = 0; f < num_factors(); ++f)
				m_lb += std::log(get_factor(f)[sub2ind(get_factor(f).vars(), m_state)]);
		}
	}

	///
	/// \brief Run the Gibbs sampler.
	///
	void run() {

		double score = 0.0;
		for (size_t f = 0; f < num_factors(); ++f)
			score += std::log(get_factor(f)[sub2ind(get_factor(f).vars(), m_state)]);

		for (size_t j = 0, i = 0; i < m_num_samples; ++i) {// keep nSamples evenly space samples
			size_t jNext = (size_t) ((1.0 + i) / m_num_samples * m_num_iter);	//   among nIter steps
			for (; j < jNext; ++j) {

				// Each iteration, go over all the variables:
				for (size_t v = 0; v < nvar(); ++v) {
					if (var(v).states() == 0)
						continue;             //   (make sure they're non-empty)
					const flist& factors = with_variable(var(v)); // get the factors they are involved with
					factor F(var(v), 1.0);
					for (flist::const_iterator f = factors.begin();
							f != factors.end(); ++f) {
						variable_set vs = get_factor(*f).vars();
						vs /= var(v);       // and condition on their neighbors
						F *= get_factor(*f).slice(vs, sub2ind(vs, m_state));
					}
					if (!m_no_best)
						score -= std::log(F[m_state[v]]);// remove current value
					if (m_temp != 1.0)
						m_state[v] = (F ^ m_temp).sample(); // then draw a new value
					else
						m_state[v] = F.sample();           //   (annealed or not)
					if (!m_no_best) {
						score += std::log(F[m_state[v]]);// update score incrementally
						if (score > m_lb) {
							m_lb = score;
							m_best_config = m_state;
						}					//   and keep the best
					}
				}						///// end iterating over each variable

				// After each sweep, update our statistics
				if (!m_no_best && isinf(score)) {// if we're still looking for a valid config
					score = 0.0; 	//   we need to update the score completely
					for (size_t f = 0; f < num_factors(); ++f)
						score += std::log(
								get_factor(f)[sub2ind(get_factor(f).vars(), m_state)]);
				}
				if (!m_no_beliefs) {		// if we're keeping marginal estimates
					for (size_t f = 0; f < num_factors(); ++f)	//   then run through the factors
						m_beliefs[f][sub2ind(get_factor(f).vars(), m_state)] += 1.0
								/ m_num_samples;	//   and update their status
				}
				if (m_temp_min != m_temp_max)
					m_temp += (m_temp_max - m_temp_min) / m_num_iter;// update temperature if annealed

			}
			m_samples.push_back(m_state);
		}
	}

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
		return belief(with_variable(v)[0]);
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
			const std::map<size_t, size_t>& old2new, const graphical_model& orig ) {

		// Open the output file
		std::ofstream out(file_name);
		if (out.fail()) {
			throw std::runtime_error("Error while opening the output file.");
		}

		// Ouput marginals (ie, beliefs)
		if (!m_no_beliefs) {
			out << "MAR" << std::endl;
			out << orig.nvar();
			for (vindex i = 0; i < orig.nvar(); ++i) {
				variable v = orig.var(i);
				try { // evidence variable
					size_t val = evidence.at(i);
					out << " " << v.states();
					for (size_t k = 0; k < v.states(); ++k) {
						out << " " << std::fixed << std::setprecision(6)
							<< (k == val ? 1.0 : 0.0);
					}
				} catch(std::out_of_range& e) { // non-evidence variable
					vindex vx = old2new.at(i);
					variable VX = var(vx);
					out << " " << VX.states();
					for (size_t j = 0; j < VX.states(); ++j)
						out << " " << std::fixed << std::setprecision(6) << belief(VX)[j];
				}
			} // end for
			out << std::endl;
		}

		// Output MAP configuration
		if (!m_no_best) {
			out << "MAP" << std::endl;
			out << orig.nvar();
			for (vindex i = 0; i < orig.nvar(); ++i) {
				try { // evidence variable
					size_t val = evidence.at(i);
					out << " " << val;
				} catch(std::out_of_range& e) { // non-evidence variable
					vindex j = old2new.at(i);
					out << " " << m_best_config[j];
				}
			}
			out << std::endl;
		}

		// Close the output file
		out.close();
	}

private:
	// Members:

	bool m_no_best;
	bool m_no_beliefs;
	size_t m_num_samples;
	size_t m_num_iter;
	std::vector<index> m_state;
	std::vector<index> m_init_state;
	variable_order_t m_order;
	std::vector<index> m_best_config;
	double m_lb;
	std::vector<factor> m_beliefs;
	std::vector<std::vector<index> > m_samples;
	double m_temp_min;
	double m_temp_max;
	double m_temp;

};

} // namespace

#endif // re-include

