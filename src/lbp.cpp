/*
 * lbp.cpp
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

/// \file lbp.cpp
/// \brief Loopy Belief Propagation (LBP) algorithm
/// \author Radu Marinescu radu.marinescu@ie.ibm.com
///

#include "lbp.h"

namespace merlin {

void lbp::init() {
	// Start the timer and store it
	m_start_time = timeSystem();

	// Prologue
	std::cout << "[LBP] + inference task   : " << "MAR" << std::endl;
	std::cout << "[LBP] + schedule         : " << m_sched << std::endl;
	std::cout << "[LBP] + exact inference  : " << "No" << std::endl;

	std::cout << "[LBP] Created factor graph with " << num_nodes()
			<< " and " << num_edges() << " edges" << std::endl;

	if (m_debug) {
		std::cout << "Variable to node (local factor) map:" << std::endl;
		for (size_t i = 0; i < m_vindex.size(); ++i) {
			std::cout << "  var " << i << " : " << m_vindex[i] << " " << get_factor(m_vindex[i]) << std::endl;
		}
		std::cout << "All node in the factor graph:" << std::endl;
		for (size_t n = 0; n < num_nodes(); ++n) {
			if (is_var_node(n)) {
				std::cout << "  node " << n << " is variable " << get_factor(n) << std::endl;
			} else {
				std::cout << "  node " << n << " is factor " << get_factor(n) << std::endl;
			}
		}
		std::cout << "Factor graph adjacencies:" << std::endl;
		for (size_t n = 0; n < num_nodes(); ++n) {
			const set<edge_id>& nbrs = neighbors(n);
			std::cout << "  node " << n << " : ";
			for (set<edge_id>::const_iterator j = nbrs.begin();
					j != nbrs.end(); ++j) {
				std::cout << *j << " ";
			}
			std::cout << std::endl;
		}
	}

	// Copy initial beliefs from factors
	m_beliefs = std::vector<factor>(m_factors);
	m_msg = std::vector<factor>();
	m_msg.resize(2 * num_edges());
	// Initialize messages to the identity
	for (size_t e = 0; e < 2 * num_edges(); ++e) {
		if (edge(e) != edge_id::NO_EDGE) {  // f'n of the right variables
			m_msg[e] = factor(get_factor(edge(e).first).vars()
							& get_factor(edge(e).second).vars(), 1.0);
		}
	}
	// Copy that as "updated" message list
	m_msg_new = std::vector<factor>(m_msg);

	m_logz = 0.0; // compute initial partition f'n estimate
	for (size_t f = 0; f < num_factors(); ++f) {
		bel(f) /= bel(f).sum(); // normalize the beliefs
		m_logz += (bel(f) * log(get_factor(f))).sum() + obj_entropy(f); // and compute the free energy estimate
	}

    // For priority scheduling
	if (m_sched == Schedule::Priority) {
		for (size_t e = 0; e < 2 * num_edges(); ++e) { //  initialize all edges to infinity
			if (edge(e) != edge_id::NO_EDGE)
				m_priority.insert(std::numeric_limits<double>::infinity(), e);
		}
	} else {
		for (size_t f = 0; f < num_factors(); ++f) {
			m_forder.push_back(f); // for fixed scheduling, get a default order
		}
	}

	if (m_debug) {
		std::cout << "Initial log partition is " << m_logz << std::endl;
		std::cout << "Initial (normalized) beliefs:" << std::endl;
		for (size_t i = 0; i < m_beliefs.size(); ++i) {
			std::cout << "  " << m_beliefs[i] << std::endl;
		}
	}
}


// Write the solution to the output stream
void lbp::write_solution(std::ostream& out, const std::map<size_t, size_t>& evidence,
		const std::map<size_t, size_t>& old2new, const graphical_model& orig,
		const std::set<size_t>& dummies, int output_format) {

	if (output_format == MERLIN_OUTPUT_JSON) {

		double val = m_logz + std::log(orig.get_global_const());

		out << "{";
		out << " \"algorithm\" : \"lbp\", ";
		out << " \"iterations\" : " << m_stop_iter * num_factors() << ", ";
		out << " \"task\" : \"MAR\", ";
		out << " \"value\" : " << std::fixed
			<< std::setprecision(MERLIN_PRECISION)
			<< (m_logz + std::log(orig.get_global_const())) << ", ";

		if (isnan(val) || isinf(val)) { // probability of evidence is 0
			out << " \"status\" : \"false\", ";
			out << " \"message\" : \"Inconsistent evidence or underflow\", ";
			out << " \"marginals\" : [] ";
		} else {
			out << " \"status\" : \"true\", ";
			out << " \"message\" : \"Consistent evidence\", ";
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
		}

		out << "}";
	} else if (output_format == MERLIN_OUTPUT_UAI) {

		double val = m_logz + std::log(orig.get_global_const());

		out << "PR" << std::endl;
		out << std::fixed << std::setprecision(MERLIN_PRECISION)
			<< m_logz + std::log(orig.get_global_const()) << " (" << std::scientific
			<< std::setprecision(MERLIN_PRECISION)
			<< std::exp(m_logz + std::log(orig.get_global_const())) << ")" << std::endl;

		if (isnan(val) || isinf(val)) {
			out << "STATUS" << std::endl;
			out << "false: Inconsistent evidence or underflow" << std::endl;
		} else {
			out << "STATUS" << std::endl;
			out << "true: Consistent evidence" << std::endl;
		}

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
				for (size_t j = 0; j < VX.states(); ++j) {
					out << " " << std::fixed
						<< std::setprecision(MERLIN_PRECISION)
						<< belief(VX)[j];
				}
			}
		} // end for
		out << std::endl;
	} else {
		std::string err_msg("Unknown output format.");
		throw std::runtime_error(err_msg);
	}
}

void lbp::run() {

	init();

	// it's easier to count updates than "iterations"
	size_t stopIter = m_stop_iter * num_factors();

	bool ok = true;
	double dObj = m_stop_obj + 1.0, dMsg = m_stop_msg + 1.0;// initialize termination values
	size_t iter = 0, print = 1; // count updates and "iterations" for printing
	size_t f, n = 0;

	std::cout << "[LBP] Begin message passing over factor graph ..." << std::endl;
	for (; dMsg >= m_stop_msg && iter < stopIter && dObj >= m_stop_obj;) {

		if (m_sched == Schedule::Priority) { // priority schedule =>
			f = edge(m_priority.top().second).second; // get next factor for update from queue
			m_priority.pop();
		} else { // fixed schedule =>
			f = m_forder[n]; // get next factor from list
			if (++n == m_forder.size()) {
				n = 0; // loop over to the beginning of the fixed order
			}
		}

		if (m_sched != Schedule::Flood) {// For non-"flood" schedules,
			factor logF = log(get_factor(f));
			dObj = 0.0; // compute new belief and update objective:
			dObj -= (belief(f) * logF).sum() + obj_entropy(f); //   remove old contribution
			accept_incoming(f); //   accept all messages into factor f
			m_logz += dObj += (belief(f) * logF).sum() + obj_entropy(f); //   re-add new contribution
		}
		update_outgoing(f);		//   update outgoing messages from factor f

		if (m_sched == Schedule::Priority) {
			dMsg = m_priority.top().first; // priority schedule => easy to check msg stop
		} else if (m_stop_msg > 0 && n == 0) { // else check once each time through all factors
			dMsg = 0.0;
			for (size_t e = 0; e < 2 * num_edges(); ++e) {
				dMsg = std::max(dMsg, m_msg_new[e].distance(m_msg[e], m_dist));
			}
		}

		if (m_sched == Schedule::Flood && n == 0) { // for flooding schedules, recalculate all
			dObj = m_logz;
			m_logz = 0.0; //   the beliefs and objective now
			for (size_t f = 0; f < num_factors(); ++f) {
				accept_incoming(f);
				m_logz += (belief(f) * log(get_factor(f))).sum() + obj_entropy(f);
			}
			dObj -= m_logz;
		}

		// check if NaN or -inf
		if (isnan(m_logz) || isinf(m_logz)) {
			ok = false;
			break;
		}

		if (iter > print * num_factors()) {
			print++;
			std::cout << "  logZ: " << std::fixed << std::setw(12)
				<< std::setprecision(MERLIN_PRECISION)
				<< m_logz << " (" << std::scientific
				<< std::setprecision(MERLIN_PRECISION)
				<< std::exp(m_logz) << ") ";
			std::cout << "\td=" << std::scientific
				<< std::setprecision(MERLIN_PRECISION)
				<< dObj << "\tm=" << dMsg << "\t time="  << std::fixed
				<< std::setprecision(MERLIN_PRECISION)
				<< (timeSystem() - m_start_time)
				<< "\ti=" << iter << std::endl;
		}

		iter++;
	}

	// Output solution (UAI output format)
	std::cout << "[LBP] Converged after " << iter << " iterations in "
		<< (timeSystem() - m_start_time) << " seconds" << std::endl;
	std::cout << "PR" << std::endl;
	std::cout << std::fixed << std::setprecision(MERLIN_PRECISION)
		<< m_logz << " (" << std::scientific
		<< std::setprecision(MERLIN_PRECISION)
		<< std::exp(m_logz) << ")" << std::endl;
	if (!ok) {
		std::cout << "STATUS" << std::endl;
		std::cout << "false: Inconsistent evidence or underflow" << std::endl;
	} else {
		std::cout << "STATUS" << std::endl;
		std::cout << "true: Consistent evidence" << std::endl;
	}
	std::cout << "MAR" << std::endl;
	std::cout << nvar();
	for (size_t v = 0; v < m_vindex.size(); ++v) {
		variable VX = var(v);
		std::cout << " " << VX.states();
		for (size_t j = 0; j < VX.states(); ++j) {
			std::cout << " " << std::fixed
				<< std::setprecision(MERLIN_PRECISION) << belief(VX)[j];
		}
	}
	std::cout << std::endl;

	if (m_debug) {
		std::cout << "Final log partition function is " << m_logz << std::endl;
		std::cout << "Final (normalized) beliefs\n";
		for (size_t i = 0; i < m_beliefs.size(); ++i) {
			std::cout << m_beliefs[i] << std::endl;
		}
	}
}

} // end namespace


