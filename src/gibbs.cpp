/*
 * gibbs.cpp
 *
 *  Created on: 12 Sep 2018
 *      Author: radu
 *
 * Copyright (c) 2018, International Business Machines Corporation
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

/// \file gibbs.cpp
/// \brief Gibbs sampling
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#include "gibbs.h"

namespace merlin {

void gibbs::init() {

	// Start the timer and store it
	m_start_time = timeSystem();
	rand_seed(); // initialize the random number generator

	// Prologue
	std::cout << "[GIBBS] + inference task   : " << m_task << std::endl;
	std::cout << "[GIBBS] + exact inference  : " << "No" << std::endl;

	// Initialize sampler with a random state
	m_samples.clear();
	m_state.resize(nvar());// might want to search for a good initialization
	for (size_t i = 0; i < nvar(); ++i) {
		m_state[i] = randi2(var(i).states());
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

void gibbs::run() {

	init();

	std::cout << "[GIBBS] Initial score: " << m_lb << std::endl;
	std::cout << "[GIBBS] Initial state: ";
	std::copy(m_state.begin(), m_state.end(),
			std::ostream_iterator<size_t>(std::cout, " "));
	std::cout << std::endl;
	std::cout << "[GIBBS] Start Gibbs sampling ... " << std::endl;

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
	std::cout << "[GIBBS] Finished in " << (timeSystem() - m_start_time) << " seconds" << std::endl;
	std::cout << "[GIBBS] Final score: " << m_lb << std::endl;
	for (size_t v = 0; v < nvar(); ++v) { //   then run through the factors
		for (size_t s = 0; s < m_samples.size(); ++s) {
			std::vector<index>& state = m_samples[s];
			m_beliefs[v][state[v]] += 1.0 / m_num_samples;
		}
		//m_beliefs[v] /= m_num_samples;
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
						<< std::setprecision(MERLIN_PRECISION) << belief(i)[j];
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
	default:
		break;
	}

}

void gibbs::write_solution(const char* file_name, const std::map<size_t, size_t>& evidence,
		const std::map<size_t, size_t>& old2new, const graphical_model& orig ) {

	// Open the output file (for overwrite or append)
	std::ofstream out(file_name);
	if (out.fail()) {
		std::string err_msg("Cannot open the output file: ");
		err_msg += std::string(file_name);
		throw std::runtime_error(err_msg);
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
						out << " " << std::fixed
							<< std::setprecision(MERLIN_PRECISION)
							<< (k == val ? 1.0 : 0.0);
					}
				} catch(std::out_of_range& e) { // non-evidence variable
					vindex vx = old2new.at(i);
					variable VX = var(vx);
					out << " " << VX.states();
					for (size_t k = 0; k < VX.states(); ++k) {
						out << " " << std::fixed
							<< std::setprecision(MERLIN_PRECISION)
							<< belief(VX)[k];
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
	default:
		break;
	}

	// Close the output file
	out.close();
}

void gibbs::write_solution(std::ostream& out, const std::map<size_t, size_t>& evidence,
		const std::map<size_t, size_t>& old2new, const graphical_model& orig,
		const std::set<size_t>& dummies, int output_format) {

	if (output_format == MERLIN_OUTPUT_JSON) {
		out << "{";
		out << " \"algorithm\" : \"gibbs\", ";
		out << " \"samples\" : " << m_samples.size() << ", ";
		switch (m_task) {
		case Task::MAR:
			{
				out << " \"task\" : \"MAR\", ";
				out << " \"value\" : " << std::fixed
					<< std::setprecision(MERLIN_PRECISION)
					<< (m_lb + std::log(orig.get_global_const())) << ", ";
				out << " \"status\" : \"true\", ";
				out << " \"marginals\" : [ ";
				for (vindex i = 0; i < orig.nvar(); ++i) {
					if (dummies.find(i) != dummies.end()) {
						continue; // skip dummy variables from output
					}

					variable v = orig.var(i);
					out << "{";
					out << " \"variable\" : " << v.label() << ", ";
					out << " \"states\" : " << v.states() << ", ";
					out << " \"probabilities\" : [";
					try { // evidence variable
						size_t val = evidence.at(i);
						for (size_t k = 0; k < v.states(); ++k) {
							out << std::fixed
								<< std::setprecision(MERLIN_PRECISION)
								<< (k == val ? 1.0 : 0.0);
							if (k != v.states() - 1) {
								out << ", ";
							}
						}
						out << "] ";
					} catch(std::out_of_range& e) { // non-evidence variable
						vindex vx = old2new.at(i);
						variable VX = var(vx);
						for (size_t k = 0; k < VX.states(); ++k) {
							out << std::fixed
								<< std::setprecision(MERLIN_PRECISION)
								<< belief(VX)[k];
							if (k != VX.states() - 1) {
								out << ", ";
							}
						}
						out << "] ";
					}

					out << "}";
					if (i != orig.nvar() - 1) {
						out << ", ";
					}
				} // end for
				out << "] ";
				break;
			}
		case Task::MAP:
			{
				out << " \"task\" : \"MAP\", ";
				out << " \"value\" : " << std::fixed
					<< std::setprecision(MERLIN_PRECISION)
					<< (m_lb + std::log(orig.get_global_const())) << ", ";
				out << " \"status\" : \"true\", ";
				out << " \"solution\" : [ ";
				for (vindex i = 0; i < orig.nvar(); ++i) {
					if (dummies.find(i) != dummies.end()) {
						continue; // skip dummy variables from output
					}

					out << "{";
					out << " \"variable\" : " << i << ",";

					try { // evidence variable
						size_t val = evidence.at(i);
						out << " \"value\" : " << val;
					} catch(std::out_of_range& e) { // non-evidence variable
						vindex j = old2new.at(i);
						out << " \"value\" : " << m_best_config[j];
					}
					out << "}";
					if (i != orig.nvar() - 1) {
						out << ", ";
					}
				}
				out << "] ";
				break;
			}
		default:
			break;
		} // end switch
	} else if (output_format == MERLIN_OUTPUT_UAI) {
		// Ouput marginals (ie, beliefs)
		switch (m_task) {
		case Task::MAR:
			{
				out << "MAR" << std::endl;
				out << orig.nvar() - dummies.size();
				for (vindex i = 0; i < orig.nvar(); ++i) {
					if (dummies.find(i) != dummies.end()) {
						continue; // skip dummy variables from output
					}

					variable v = orig.var(i);
					try { // evidence variable
						size_t val = evidence.at(i);
						out << " " << v.states();
						for (size_t k = 0; k < v.states(); ++k) {
							out << " " << std::fixed
								<< std::setprecision(MERLIN_PRECISION)
								<< (k == val ? 1.0 : 0.0);
						}
					} catch(std::out_of_range& e) { // non-evidence variable
						vindex vx = old2new.at(i);
						variable VX = var(vx);
						out << " " << VX.states();
						for (size_t k = 0; k < VX.states(); ++k) {
							out << " " << std::fixed
								<< std::setprecision(MERLIN_PRECISION)
								<< belief(VX)[k];
						}
					}
				} // end for
				out << std::endl;
				break;
			}
		case Task::MAP:
			{
				out << "MAP" << std::endl;
				out << orig.nvar() - dummies.size();
				for (vindex i = 0; i < orig.nvar(); ++i) {
					if (dummies.find(i) != dummies.end()) {
						continue; // skip dummy variables from output
					}

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
		default:
			break;
		} // end switch
	} else {
		std::string err_msg("Unknown output format");
		throw std::runtime_error(err_msg);
	}
}


} // end namespace



