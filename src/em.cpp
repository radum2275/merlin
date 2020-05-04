/*
 * em.cpp
 *
 *  Created on: Aug 2, 2013
 *      Author: radu
 *
 * Copyright (c) 2019, International Business Machines Corporation
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

/// \file em.cpp
/// \brief EM parameter learning for Bayes nets (directed models)
/// \author Radu Marinescu radu.marinescu@ie.ibm.com
///

#include "em.h"
#include "cte.h"

namespace merlin {

// The E step (expectation)
void em::e_step() {

	if (m_debug) {
		std::cout << "[DEBUG] Begin E-step" << std::endl;
		std::cout << "[DEBUG] Reseting counts" << std::endl;
	}

	// Reset the counts to 0
	for (size_t i = 0; i < m_counts.size(); ++i) {
		m_counts[i].fill(0.0);
	}

	// Reset the log-likelihood
	m_loglikelihood = 0.0;

	// Go over each of the training examples and compute (expected) counts
	size_t M = m_dataset.size();
	for (size_t m = 0; m < M; ++m) {
		std::vector<observation>& om = m_dataset[m];
		std::vector<int> evidence(om.size(), -1);
		std::vector<observation> virtualEvidence;
		std::set<vindex> dummies;
		for (size_t i = 0; i < om.size(); ++i) {
			observation obs = om[i];
			if (obs.is_virtual()) {
				virtualEvidence.push_back(obs);
			} else if (obs.is_observed()) {
				evidence[i] = obs.val();
			}
		}

		// If virtual evidence is present, extend the model with extra variables
		if (!virtualEvidence.empty()) {
			graphical_model gm(m_gmo); // copy the current graphical model
			vindex idx = gm.nvar();	// last variable index
			std::vector<factor> nfs = m_gmo.get_factors();
			for (size_t j = 0; j < virtualEvidence.size(); ++j) {
				observation obs = virtualEvidence[j];
				vindex x = (vindex)obs.var();		// virtual evidence variable
				likelihood l = obs.likelihood();	// likelihood
				variable xvar = gm.var(x);
				variable uvar(idx++, 2);
				variable_set vs;
				vs |= xvar;
				vs |= uvar;
				factor f(vs, 0.0);
				f.set_child(uvar.label());
				evidence.push_back(0); // first value of the U variable
				dummies.insert(uvar.label()); // keep track of dummy variables
				for (size_t k = 0; k < l.size(); ++k) {
					f.set(k, l[k]);
					f.set(k + xvar.states(), 1.0 - l[k]);
				}

				gm.add_factor(f); // add the extra factor
			}

			// Create a new instance of the inference engine with the extra factors
			cte temp(gm);
			temp.set_properties(m_properties);
			temp.init();
			bool ok = temp.propagate_evidence(evidence);
			if (ok) { // silently ignore inconsistent evidence!
				m_loglikelihood += temp.logZ();
				for (size_t i = 0; i < m_counts.size(); ++i) {
					variable_set vs = m_counts[i].vars();
					temp.joint_marginal(vs, evidence);
					m_counts[i] += temp.get_joint_marginal();
				}
			}
		} else { // regular evidence

			// Propagate the regular evidence
			bool ok = m_infer.propagate_evidence(evidence);
			if (ok) { // silently ignore inconsistent evidence!
				m_loglikelihood += m_infer.logZ(); // probability of evidence
				for (size_t i = 0; i < m_counts.size(); ++i) {
					variable_set vs = m_counts[i].vars();
					m_infer.joint_marginal(vs, evidence);
					m_counts[i] += m_infer.get_joint_marginal();
				}
			}
		}
	}

	if (m_debug) {
		std::cout << "[DEBUG] End E-step" << std::endl;
		std::cout << "[DEBUG] Log-likelihood = " << m_loglikelihood << std::endl;
		std::cout << "[DEBUG] Family Counts:" << std::endl;
		for (size_t i = 0; i < m_counts.size(); ++i) {
			std::cout << "[DEBUG] " << m_counts[i] << std::endl;
		}
	}
}

// The M step (maximization)
void em::m_step() {

	if (m_debug) {
		std::cout << "[DEBUG] Begin M-step" << std::endl;
	}

	// Update the new parameters by dividing the corresponding counts
	const std::vector<factor>& thetas = m_gmo.get_factors();
	std::vector<factor> new_thetas(thetas.size());

	// counts factors must be aligned to the theta factors (by construction)
	for (size_t i = 0; i < thetas.size(); ++i) {
		const factor& th = thetas[i];
		int x = th.get_child();
		assert(x >= 0);

		variable_set vx(m_gmo.var(x));
		factor sum = m_counts[i].sum(vx);
		factor temp = m_counts[i] / sum;
		variable_set scope = th.vars();
		factor new_th = th; // copy the previous factor
		index_config cv1(scope, true);
		config_index cv2(sum.vars(), true);
		for (size_t j = 0; j < th.num_states(); ++j) {
			std::map<size_t, size_t> config = cv1.convert(j);
			size_t k = cv2.convert(config);
			if (sum.get(k) != 0) {
				double v = temp.get(j);
				new_th.set(j, v);
			}
		}

		new_thetas[i] = new_th;
	}

	// Set the new parameters (discard the old ones)
	m_gmo = graphical_model(new_thetas);
	m_infer.reinit(m_gmo.get_factors());
//	m_infer = cte(m_gmo);
//	m_infer.set_properties(m_properties);
//	m_infer.init(); // re-initialize the join tree

	if (m_debug) {
		std::cout << "[DEBUG] End M-step" << std::endl;
		std::cout << "[DEBUG] Updated parameters:" << std::endl;
		for (size_t i = 0; i < new_thetas.size(); ++i) {
			std::cout << "[DEBUG]  " << new_thetas[i] << std::endl;
		}
	}
}

void em::init() {

	// Prologue
	std::cout << "[EM] + inferene method  : " << "CTE" << std::endl;
	std::cout << "[EM] + iterations       : " << m_iterations << std::endl;
	std::cout << "[EM] + epsilon          : " << m_epsilon << std::endl;
	std::cout << "[EM] + ordering method  : " << m_order_method << std::endl;

	// Collect the missing values stats
	double missingVals = 0.0, virtualVals = 0.0, allVals = 0.0;
	for (size_t i = 0; i < m_dataset.size(); ++i) {
		for (size_t j = 0; j < m_dataset[i].size(); ++j) {
			allVals += 1;
			if (!m_dataset[i][j].is_observed() &&
				!m_dataset[i][j].is_virtual()) {

				missingVals += 1;
			} else if (m_dataset[i][j].is_virtual()) {
				virtualVals += 1;
			}
		}
	}

	// Percentage of missing values
	double perc1 = (missingVals/allVals)*100.0;
	double perc2 = (virtualVals/allVals)*100.0;
	std::cout << "[EM] + dataset size     : " << m_dataset.size() << " examples" << std::endl;
	std::cout << "[EM] + missing values   : " << missingVals << "/" << allVals << " (" << perc1 << "%)" << std::endl;
	std::cout << "[EM] + virtual evidence : " << virtualVals << "/" << allVals << " (" << perc2 << "%)" << std::endl;

	// Set the CTE properties
	std::ostringstream oss;
	oss << "Order=MinFill" << ","
		<< "Task=MAR" << ","
		<< "Debug=" << (m_debug ? "1" : "0") << ","
		<< "Verbose=0";
	m_properties = oss.str();

	// Initialize the CPTs uniformly at random
	if (m_init_method == InitMethod::Uniform) {
		m_gmo.uniform_bayes();
	} else if (m_init_method == InitMethod::Random){
		m_gmo.random_bayes();
	}

	// Initialize the join tree
	m_infer = cte(m_gmo);
	m_infer.set_properties(m_properties);
	m_infer.init(); // initialize the join tree

	// Initialize the families (Bayes nets only)
	size_t n = m_gmo.nvar();
	const std::vector<factor>& factors = m_gmo.get_factors();
	size_t m = factors.size();
	m_families.resize(n);
	m_counts.resize(m);
	for (size_t i = 0; i < m; ++i) {
		const factor& f = factors[i];
		int child = f.get_child();
		assert(child >= 0); // Bayes factor (CPT)
		variable_set ps; // the parents set
		std::vector<vindex> pa;
		for (variable_set::const_iterator it = f.vars().begin();
				it != f.vars().end(); ++it) {
			if (it->label() != (size_t)child) {
				pa.push_back(it->label());
				ps |= *it;
			}
		}

		m_families[child] = pa;
		m_counts[i] = factor(f.vars(), 0.0);
	}

	if (m_debug) {
		std::cout << "Families:" << std::endl;
		for (size_t i = 0; i < m_families.size(); ++i) {
			std::cout << "var " << i << ": ";
			std::copy(m_families[i].begin(), m_families[i].end(),
					std::ostream_iterator<int>(std::cout, " "));
			std::cout << std::endl;
		}

		std::cout << "Initial parameters:" << std::endl;
		for (size_t i = 0; i < m_gmo.get_factors().size(); ++i) {
			std::cout << " " << m_gmo.get_factor(i) << std::endl;
		}
	}

	std::cout << std::endl;
}

// The EM algorithm
void em::run() {

	// Initialize the algorithm
	init();
	std::cout << "[EM] Begin parameter learning ..." << std::endl;
	double timeStart = timeSystem();

	// Iterate the e-step and the m-step over the dataset
	double loglikelihood = 0;
	for (int i = 0; i < m_iterations; ++i) {
		e_step();
		m_step();
		std::cout << " " << i << ": log-likelihood = " << m_loglikelihood << std::endl;

		double delta = std::fabs(m_loglikelihood - loglikelihood);
		if (delta <= m_epsilon) {
			std::cout << "[EM] Converged to log-likelihood "
					<< m_loglikelihood << " after " << i << " iterations" << std::endl;
			break; // converged skip this for now
		}

		loglikelihood = m_loglikelihood;
	}

	std::cout << "[EM] Finished parameter learning" << std::endl;
	std::cout << "[EM] Time elapsed is " << (timeSystem() - timeStart) << " seconds" << std::endl;
}

// Write the model file with the new parameters learned from the training dataset
void em::write_solution(std::ostream& out, const graphical_model& orig) {

	if (orig.is_bayes()) {
		m_gmo.write_bayes(out);
	} else {
		m_gmo.write(out);
	}
}


} // end namespace

