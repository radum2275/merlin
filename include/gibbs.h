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
	gibbs* clone() const {
		gibbs* fg = new gibbs(*this);
		return fg;
	}

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
	void init() {

		// Start the timer and store it
		m_start_time = timeSystem();

		// Prologue
		std::cout << VERSIONINFO << std::endl << COPYRIGHT << std::endl;
		std::cout << "Initialize inference engine ..." << std::endl;
		std::cout << "+ tasks supported  : PR, MAR, MAP" << std::endl;
		std::cout << "+ algorithm        : " << "IJGP" << std::endl;
		std::cout << "+ inference task   : " << m_task << std::endl;

		// Initialize sampler with a random state
		m_samples.clear();
		m_state.resize(nvar());// might want to search for a good initialization
		for (size_t i = 0; i < nvar(); ++i) {
			m_state[i] = randi(var(i).states());
		}

		// Initialize beliefs
		m_beliefs.resize(nvar());
		for (size_t i = 0; i < nvar(); ++i) {
			m_beliefs[i] = factor(var(i), 0.0);
		}

		// Initialize best configuration
		m_best_config = m_state;
		m_lb = 0.0;
		for (size_t f = 0; f < num_factors(); ++f) {
			m_lb += std::log(get_factor(f)[sub2ind(get_factor(f).vars(), m_state)]);
		}
	}

	///
	/// \brief Run the Gibbs sampler.
	///
	void run() {

		init();

		std::cout << "Initial score: " << m_lb << std::endl;
		std::cout << "Initial state: ";
		std::copy(m_state.begin(), m_state.end(),
				std::ostream_iterator<size_t>(std::cout, " "));
		std::cout << std::endl;

		// Generate the samples
		double score = m_lb;
		for (size_t j = 0, i = 0; i < m_num_samples; ++i) {// keep nSamples evenly space samples
			size_t jNext = (size_t) ((1.0 + i) / m_num_samples * m_num_iter);	//   among nIter steps
			for (; j < jNext; ++j) {

				// Each iteration, go over all the variables:
				std::vector<index> sample(nvar(), -1); // new sample to be generated
				for (size_t v = 0; v < nvar(); ++v) {
					assert(var(v).states() != 0); // (make sure they're non-empty)

					const flist& factors = with_variable(var(v)); // get the factors they are involved with
					factor F(var(v), 1.0);
					for (flist::const_iterator f = factors.begin();
							f != factors.end(); ++f) {
						variable_set vs = get_factor(*f).vars();
						vs /= var(v);       // and condition on their neighbors
						F *= get_factor(*f).slice(vs, sub2ind(vs, m_state));
					}

					score -= std::log(F[m_state[v]]);// remove current value
					if (m_temp != 1.0) {
						sample[v] = (F ^ m_temp).sample(); // then draw a new value
					} else {
						sample[v] = F.sample();           //   (annealed or not)
					}

					score += std::log(F[sample[v]]);// update score incrementally
				} // end iterating over each variable

				if (score > m_lb) {
					m_lb = score;
					m_best_config = sample;
				} //   and keep the best

				// We have a new sample
				m_state = sample;

				// After each sweep, update our statistics
				if (isinf(score)) {// if we're still looking for a valid config
					score = 0.0; 	//   we need to update the score completely
					for (size_t f = 0; f < num_factors(); ++f)
						score += std::log(
								get_factor(f)[sub2ind(get_factor(f).vars(), m_state)]);
				}
//				for (size_t v = 0; v < nvar(); ++v) { //   then run through the factors
//					m_beliefs[v][m_state[v]] += 1.0 / m_num_samples; // and update their status
//				}

				if (m_temp_min != m_temp_max)
					m_temp += (m_temp_max - m_temp_min) / m_num_iter;// update temperature if annealed

			}

			m_samples.push_back(m_state);
		}

		// check out the samples
		if (m_debug) {
			std::cout << "Samples generated: " << m_samples.size() << std::endl;
			for (size_t s = 0; s < m_samples.size(); ++s) {
				std::copy(m_samples[s].begin(), m_samples[s].end(),
					std::ostream_iterator<index>(std::cout, " "));
				std::cout << std::endl;
			}
		}

		// update the beliefs
		for (size_t v = 0; v < nvar(); ++v) { //   then run through the factors
			for (size_t s = 0; s < m_samples.size(); ++s) {
				std::vector<index>& state = m_samples[s];
				m_beliefs[v][state[v]] += 1.0;
			}
			m_beliefs[v] /= m_num_samples;
		}


		// Ouput marginals (ie, beliefs)
		switch (m_task) {
		case Task::MAR:
			{
				std::cout << "MAR" << std::endl;
				std::cout << nvar();
				for (vindex i = 0; i < nvar(); ++i) {
					variable VX = var(i);
					std::cout << " " << VX.states();
					for (size_t j = 0; j < VX.states(); ++j) {
						std::cout << " " << std::fixed
							<< std::setprecision(6) << belief(i)[j];
					}
				} // end for
				std::cout << std::endl;
				break;
			}
		case Task::MAP:
			{
				std::cout << "MAP" << std::endl;
				std::cout << nvar();
				for (vindex i = 0; i < nvar(); ++i) {
					std::cout << " " << m_best_config[i];
				}
				std::cout << std::endl;
				break;
			}
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
		switch (m_task) {
		case Task::MAR:
			{
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
						for (size_t j = 0; j < VX.states(); ++j) {
							out << " " << std::fixed << std::setprecision(6) << belief(i)[j];
						}
					}
				} // end for
				out << std::endl;
				break;
			}
		case Task::MAP:
			{
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
				break;
			}
		}

		// Close the output file
		out.close();
	}

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

