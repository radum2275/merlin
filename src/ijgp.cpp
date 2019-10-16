/*
 * ijgp.cpp
 *
 *  Created on: 12 Sep 2018
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

/// \file ijgp.cpp
/// \brief Iterative Join Graph Propagation (IJGP) algorithm
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#include "ijgp.h"

namespace merlin {

// Write the solution to the output stream
void ijgp::write_solution(std::ostream& out, const std::map<size_t, size_t>& evidence,
		const std::map<size_t, size_t>& old2new, const graphical_model& orig,
		const std::set<size_t>& dummies, int output_format) {

	if (output_format == MERLIN_OUTPUT_JSON) {
		out << "{";
		out << " \"algorithm\" : \"ijgp\", ";
		out << " \"ibound\" : " << m_ibound << ", ";
		out << " \"iterations\" : " << num_iter << ", ";
		switch (m_task) {
		case Task::MAR:
			{
				double val = m_logz + std::log(orig.get_global_const());
				double prob = std::exp(val);

				out << " \"task\" : \"MAR\", ";
				out << " \"value\" : " << std::fixed
					<< std::setprecision(MERLIN_PRECISION)
					<< (m_logz + std::log(orig.get_global_const())) << ", ";
				if (prob == 0.0) { // probability of evidence is 0
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

				break;
			}
		case Task::MAP:
			{
				out << " \"task\" : \"MAP\", ";
				out << " \"value\" : " << std::fixed
					<< std::setprecision(MERLIN_PRECISION)
					<< (m_logz + std::log(orig.get_global_const())) << ", ";
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
						out << " \"value\" : " << val << ",";
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
		}
		out << "}";
	} else if (output_format == MERLIN_OUTPUT_UAI) {
		switch (m_task) {
		case Task::MAR:
			{
				double val = m_logz + std::log(orig.get_global_const());
				double prob = std::exp(val);

				out << "PR" << std::endl;
				out << std::fixed << std::setprecision(MERLIN_PRECISION)
					<< m_logz << " (" << std::scientific
					<< std::setprecision(MERLIN_PRECISION)
					<< std::exp(m_logz) << ")" << std::endl;

				if (prob == 0.0) {
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
		}
	} else {
		std::string err_msg("Unknown output format");
		throw std::runtime_error(err_msg);
	}
}

void ijgp::run() {

	// Start the timer and store it
	m_start_time = timeSystem();

	init();
	propagate(num_iter);

	// Output solution (UAI output format)
	std::cout << "[IJGP] Converged after " << num_iter << " iterations in "
			<< (timeSystem() - m_start_time) << " seconds" << std::endl;

	switch (m_task) {
	case Task::PR:
	case Task::MAR:
		{
			std::cout << "PR" << std::endl;
			std::cout << std::fixed
				<< std::setprecision(MERLIN_PRECISION)
				<< m_logz << " (" << std::scientific
				<< std::setprecision(MERLIN_PRECISION)
				<< std::exp(m_logz) << ")" << std::endl;

			if (isinf(m_logz)) {
				std::cout << "STATUS" << std::endl;
				std::cout << "false: Inconsistent evidence or underflow" << std::endl;
			} else {
				std::cout << "STATUS" << std::endl;
				std::cout << "true: Consistent evidence" << std::endl;
			}

			std::cout << "MAR" << std::endl;
			std::cout << m_gmo.nvar();
			for (vindex v = 0; v < m_gmo.nvar(); ++v) {
				variable VX = m_gmo.var(v);
				std::cout << " " << VX.states();
				for (size_t j = 0; j < VX.states(); ++j) {
					std::cout << " " << std::fixed
						<< std::setprecision(MERLIN_PRECISION)
						<< belief(VX)[j];
				}
			}
			std::cout << std::endl;

			break;
		}
	case Task::MAP:
		{
			m_lb = m_gmo.logP(m_best_config);
			std::cout << "Final Lower Bound is " << std::fixed
				<< std::setw(12) << std::setprecision(MERLIN_PRECISION)
				<< m_lb << " (" << std::scientific
				<< std::setprecision(MERLIN_PRECISION)
				<< std::exp(m_lb) << ")" << std::endl;
			std::cout << "MAP" << std::endl;
			std::cout << m_gmo.nvar();
			for (vindex v = 0; v < m_gmo.nvar(); ++v) {
				std::cout << " " << m_best_config[v];
			}
			std::cout << std::endl;

			break;
		}
	}
}

void ijgp::init() {

	// Prologue
	std::cout << "[IJGP] + i-bound          : " << m_ibound << std::endl;
	std::cout << "[IJGP] + iterations       : " << num_iter << std::endl;
	std::cout << "[IJGP] + inference task   : " << m_task << std::endl;
	std::cout << "[IJGP] + ordering method  : " << m_order_method << std::endl;
	std::cout << "[IJGP] + elimination      : ";

	if (m_order.size() == 0) { // if we need to construct an elimination ordering
		m_order = m_gmo.order(m_order_method);
		m_parents.clear(); // (new elim order => need new pseudotree) !!! should do together
		std::copy(m_order.begin(), m_order.end(),	std::ostream_iterator<size_t>(std::cout, " "));
	}
	if (m_parents.size() == 0) {     // if we need to construct a pseudo-tree
		m_parents = m_gmo.pseudo_tree(m_order);
	}

	std::cout << std::endl;
	size_t wstar = m_gmo.induced_width(m_order);
	std::cout << "[IJGP] + induced width    : " << wstar << std::endl;
	std::cout << "[IJGP] + exact inference  : " << (m_ibound >= wstar ? "Yes" : "No") << std::endl;
	std::cout << "[IJGP] + ordering time    : " << (timeSystem() - m_start_time) << " seconds" << std::endl;
	if (m_ibound >= wstar) {
		num_iter = 1; // exact inference requires 1 iteration over the join-tree
	}

	// Get the factors scopes
	vector<variable_set> fin;
	for (vector<factor>::const_iterator i = m_gmo.get_factors().begin();
			i != m_gmo.get_factors().end(); ++i) {
		fin.push_back((*i).vars());
	}

	if (m_debug) {
		std::cout << "[DEBUG] Original factors scopes:" << std::endl;
		for (size_t i = 0; i < fin.size(); ++i) {
			std::cout << i << ": " << fin[i] << std::endl;
		}
	}

	// Mark factors depending on variable i
	vector<flist> vin;
	for (size_t i = 0; i < m_gmo.nvar(); ++i) {
		vin.push_back(m_gmo.with_variable(var(i)));
	}

	vector<flist> Orig(m_gmo.num_factors()); 	// origination info: which original factors are
	for (size_t i = 0; i < Orig.size(); ++i)
		Orig[i] |= i;    					// included for the first time, and which newly
	vector<flist> New(m_gmo.num_factors()); 	// created clusters feed into this cluster

	// Initialize join-graph by running mini-buckets schematically
	std::cout << "[IJGP] Initializing join-graph ... " << std::endl;

	m_clusters.resize(m_order.size());
	for (variable_order_t::const_iterator x = m_order.begin(); x != m_order.end(); ++x) {

		if (m_debug) {
			std::cout << "  - create bucket/cluster for var " << *x << std::endl;
		}

		variable VX = var(*x);
		if (*x >= vin.size() || vin[*x].size() == 0)
			continue;  // check that we have some factors over this variable

		flist ids = vin[*x];  // list of factor IDs contained in this bucket

		// Select allocation into mini-buckets
		typedef flist::const_iterator flistIt;
		typedef std::pair<double, sPair> _INS;
		std::multimap<double, sPair> scores;
		std::map<sPair, std::multimap<double, sPair>::iterator> reverseScore;

		// Populate list of pairwise scores for aggregation
		for (flistIt i = ids.begin(); i != ids.end(); ++i) {
			for (flistIt j = ids.begin(); j != i; ++j) {
				double err = score(fin, VX, *i, *j);
				sPair sp(*i, *j);
				reverseScore[sp] = scores.insert(_INS(err, sp)); // save score
			}
			reverseScore[sPair(*i, *i)] = scores.insert(
					_INS(-1, sPair(*i, *i)));       // mark self index at -1
		}

		//// Run through until no more pairs can be aggregated:
		//   Find the best pair (ii,jj) according to the scoring heuristic and join
		//   them as jj; then remove ii and re-score all pairs with jj
		for (;;) {
			std::multimap<double, sPair>::reverse_iterator top =
					scores.rbegin();
			if (top->first < 0)
				break;                         // if can't do any more, quit
			else {
				size_t ii = top->second.first, jj = top->second.second;
				//std::cout<<"Joining "<<ii<<","<<jj<<"; size "<<(fin[ii].vars()+fin[jj].vars()).nrStates()<<"\n";
				fin[jj] |= fin[ii];                        // combine into j
				erase(vin, ii, fin[ii]);
				fin[ii] = variable_set();  //   & remove i

				Orig[jj] |= Orig[ii];
				Orig[ii].clear(); // keep track of list of original factors in this cluster
				New[jj] |= New[ii];
				New[ii].clear(); //  list of new "message" clusters incoming to this cluster

				for (flistIt k = ids.begin(); k != ids.end(); ++k) { // removing entry i => remove (i,k) for all k
					scores.erase(reverseScore[sPair(ii, *k)]);
				}
				ids /= ii;

				for (flistIt k = ids.begin(); k != ids.end(); ++k) { // updated j; rescore all pairs (j,k)
					if (*k == jj)
						continue;
					double err = score(fin, VX, jj, *k);
					sPair sp(jj, *k);
					scores.erase(reverseScore[sp]);    // change score (i,j)
					reverseScore[sp] = scores.insert(_INS(err, sp));  //
				}
			}
		}

		if (m_debug) {
			std::cout << "  - mini-buckets: " << ids.size() << std::endl;
		}

		// Eliminate individually each mini-bucket
		vector<findex> alphas;
		for (flistIt i = ids.begin(); i != ids.end(); ++i) {
			//
			// Create new cluster alpha over this set of variables; save function parameters also
			findex alpha = findex(-1);
			alpha = add_factor(factor(fin[*i]));
			alphas.push_back(alpha);
			m_clusters[*x] |= alpha;
			m_cluster2var[alpha] = *x;

			//std::cout << "  " << alpha << " : " << fin[*i] << std::endl;

			fin[*i] = fin[*i] - VX;

			// add inter clusters edges
			for (flistIt j = New[*i].begin(); j != New[*i].end(); ++j) {
				add_edge(*j, alpha);
				m_schedule.push_back(std::make_pair(*j, alpha));
			}

			// keep track of original factors
			m_originals.push_back(flist());
			m_originals[alpha] |= Orig[*i];

			// now incoming nodes to *i is just alpha
			Orig[*i].clear();
			New[*i].clear();
			New[*i] |= alpha;

			// recompute and update adjacency
			insert(vin, *i, fin[*i]);
		}

		// add extra edges between the clusters corresp. to a variable
		for (size_t i = 0; i < alphas.size() - 1; ++i) {
			add_edge(alphas[i], alphas[i+1]);
			m_schedule.push_back(std::make_pair(alphas[i], alphas[i+1]));
		}
	} // end for: variable elim order

	if (m_debug) {
		std::cout << "  - final number of clique factors is: " << m_factors.size() << std::endl;
		std::cout << "Finished initializing the join-graph." << std::endl;

		std::cout << "Propagation schedule:" << std::endl;
		for (size_t i = 0; i < m_schedule.size(); ++i) {
			std::cout << "msg (" << m_schedule[i].first << ","
					<< m_schedule[i].second << ")" << std::endl;
		}
	}

	// Create separators and cluster scopes
	size_t C = m_factors.size(), max_sep_size = 0, max_clique_size = 0;
	m_separators.resize(C);
	for (size_t i = 0; i < C; ++i) m_separators[i].resize(C);
	m_scopes.resize(C);
	for (size_t i = 0; i < C; ++i) {
		m_scopes[i] = m_factors[i].vars();
		max_clique_size = std::max(max_clique_size, m_scopes[i].size());
	}
	const std::vector<edge_id>& elist = edges();
	for (size_t i = 0; i < elist.size(); ++i) {
		findex a,b;
		a = elist[i].first;
		b = elist[i].second;
		if (a > b) continue;
		variable_set sep = m_factors[a].vars() & m_factors[b].vars();
		m_separators[a][b] = sep;
		m_separators[b][a] = sep;
		max_sep_size = std::max(max_sep_size, sep.size());
	}

	// Create incoming and outgoing lists for each cluster
	m_in.resize(C);
	m_out.resize(C);
	for (vector<std::pair<findex, findex> >::const_iterator i = m_schedule.begin();
			i != m_schedule.end(); ++i) {
		findex from = (*i).first;
		findex to = (*i).second;
		m_in[to] |= from;
		m_out[from] |= to;
	}

	// Initialize the root cluster(s)
	for (size_t i = 0; i < m_out.size(); ++i) {
		if ( m_out[i].empty() )
			m_roots |= i;
	}

	// Initialize the forward and backward messages
	size_t N = m_schedule.size();
	m_forward.resize(N);
	m_backward.resize(N);
	m_edge_indeces.resize(C);
	for (size_t i = 0; i < C; ++i) m_edge_indeces[i].resize(C);
	for (size_t i = 0; i < N; ++i) {
		size_t from = m_schedule[i].first;
		size_t to = m_schedule[i].second;
		m_edge_indeces[from][to] = i;
	}

	// Initialize the clique potentials
	for (size_t i = 0; i < m_factors.size(); ++i) {
		m_factors[i] = factor(1.0); // init

		// Compute the clique potential by multiplying he originals
		for (flist::const_iterator j = m_originals[i].begin();
				j != m_originals[i].end(); ++j) {
			m_factors[i] *= m_gmo.get_factor(*j);
		}
	}

	// Initialize beliefs (marginals)
	m_logz = 0;
	m_beliefs.clear();
	m_beliefs.resize(m_gmo.nvar(), factor(1.0));
	m_best_config.resize(m_gmo.nvar(), -1);

	// Output summary of initialization
	std::cout << "[IJGP] Created join graph with " << C << " clique factors" << std::endl;
	std::cout << "[IJGP] Number of cliques  : " << C << std::endl;
	std::cout << "[IJGP] Number of edges    : " << elist.size() << std::endl;
	std::cout << "[IJGP] Max clique size    : " << max_clique_size << std::endl;
	std::cout << "[IJGP] Max separator size : " << max_sep_size << std::endl;
	std::cout << "[IJGP] Finished initialization in " << (timeSystem() - m_start_time) << " seconds" << std::endl;

	if (m_debug) {
		std::cout << "[MERLIN DEBUG]\n";
		std::cout << "[DBG] Join-graph with " << m_factors.size() << " clusters and "
				<< elist.size() << " edges" << std::endl;
		for (size_t i = 0; i < elist.size(); ++i) {
			findex a,b;
			a = elist[i].first;
			b = elist[i].second;
			if (a > b) continue;
			std::cout << "  edge from "
					<< m_scopes[a] << " to "
					<< m_scopes[b] << " (a=" << a << ", b=" << b << ")"
					<< " sep: " << m_separators[a][b]
					<< std::endl;
		}

		std::cout << "[DBG] Forward propagation schedule:" << std::endl;
		for (size_t i = 0; i < m_schedule.size(); ++i) {
			std::cout << " msg " << m_schedule[i].first << " --> "
					<< m_schedule[i].second << std::endl;
		}
		std::cout << "[DBG] Backward propagation schedule:" << std::endl;
		vector<std::pair<findex, findex> >::reverse_iterator ri = m_schedule.rbegin();
		for (; ri != m_schedule.rend(); ++ri) {
			std::cout << " msg " << ri->second << " --> "
					<< ri->first << std::endl;
		}

		std::cout << "[DBG] Original factors per cluster:" << std::endl;
		for (size_t i = 0; i < m_originals.size(); ++i) {
			std::cout << " cl " << i << " : ";
			std::copy(m_originals[i].begin(), m_originals[i].end(),
					std::ostream_iterator<int>(std::cout, " "));
			std::cout << std::endl;
		}

		// _in and _out lists
		std::cout << "[DBG] _IN list:" << std::endl;
		for (size_t i = 0; i < m_in.size(); ++i) {
			std::cout << "  _in[" << i << "] = ";
			std::copy(m_in[i].begin(), m_in[i].end(),
					std::ostream_iterator<int>(std::cout, " "));
			std::cout << std::endl;
		}
		std::cout << "[DBG] _OUT list:" << std::endl;
		for (size_t i = 0; i < m_out.size(); ++i) {
			std::cout << "  _out[" << i << "] = ";
			std::copy(m_out[i].begin(), m_out[i].end(),
					std::ostream_iterator<int>(std::cout, " "));
			std::cout << std::endl;
		}
		std::cout << "[DBG] _ROOTS: ";
		std::copy(m_roots.begin(), m_roots.end(),
				std::ostream_iterator<int>(std::cout, " "));
		std::cout << std::endl;

		// clique factors, forward and backward
		std::cout << "[DBG] clique_factors:" << std::endl;
		for (size_t i = 0; i < m_factors.size(); ++i) {
			std::cout << "[" << i << "]: " << m_factors[i] << std::endl;
		}
		std::cout << "[DBG] _forward messages (top-down):" << std::endl;
		for (size_t i = 0; i < m_forward.size(); ++i) {
			std::cout << "(" << i << "): " << m_forward[i] << std::endl;
		}
		std::cout << "[DBG] _backward messages (bottop-up):" << std::endl;
		for (size_t i = 0; i < m_backward.size(); ++i) {
			std::cout << "(" << i << "): " << m_backward[i] << std::endl;
		}
		std::cout << "[MERLIN DEBUG]\n";
	}
}

factor ijgp::calc_belief(findex a) {

	factor bel = m_factors[a];

	// forward messages to cluster 'a'
	for (flist::const_iterator ci = m_in[a].begin();
			ci != m_in[a].end(); ++ci) {
		findex p = (*ci);
		size_t j = m_edge_indeces[p][a];
		bel *= m_forward[j];
	}

	// backward message to cluster 'a'
	for (flist::const_iterator ci = m_out[a].begin();
			ci != m_out[a].end(); ++ci) {
		findex p = (*ci);
		size_t j = m_edge_indeces[a][p];
		bel *= m_backward[j];
	}

	return bel;
}

factor ijgp::calc_belief(findex a, size_t b) {

	factor bel = m_factors[a];

	// forward messages to cluster 'a'
	for (flist::const_iterator ci = m_in[a].begin();
			ci != m_in[a].end(); ++ci) {
		findex p = (*ci);
		if (p == b) continue;
		size_t j = m_edge_indeces[p][a];
		bel *= m_forward[j];
	}

	// backward message to cluster 'a'
	for (flist::const_iterator ci = m_out[a].begin();
			ci != m_out[a].end(); ++ci) {
		findex p = (*ci);
		if (p == b) continue;
		size_t j = m_edge_indeces[a][p];
		bel *= m_backward[j];
	}

	return bel;
}

factor ijgp::incoming(findex a) {

	factor bel = m_factors[a];

	// forward messages to 'a'
	for (flist::const_iterator ci = m_in[a].begin();
			ci != m_in[a].end(); ++ci) {
		findex p = (*ci);
		size_t j = m_edge_indeces[p][a];
		bel *= m_forward[j];
	}

	return bel;
}

void ijgp::forward() {

	// Compute and propagate forward messages (top-down)
	if (m_debug) {
		std::cout << "Begin forward (top-down) pass ..." << std::endl;
	}

	m_logz = 0; // reset the log partition function
	vector<std::pair<findex, findex> >::iterator fi = m_schedule.begin();
	for (; fi != m_schedule.end(); ++fi ) {

		// Compute forward message m(a->b)
		findex a = (*fi).first;  // cluster index a
		findex b = (*fi).second; // cluster index b
		size_t ei = m_edge_indeces[a][b]; // edge index

		// variables to eliminate
		variable_set VX = m_scopes[a] - m_separators[a][b];

		// compute the belief at a (excluding message b->a)
		factor bel = calc_belief(a, b);
		m_forward[ei] = elim(bel, VX);
		m_forward[ei].normalize();
//		double mx = m_forward[ei].max(); // normalize for stability
//		m_forward[ei] /= mx;
//		m_logz += std::log(mx);

		if (m_debug) {
			std::cout << " - Sending forward msg from " << a << " to " << b;
			std::cout << "  - forward msg (" << a << "," << b << "): elim = " << VX << std::endl;
			std::cout << "  -> " << m_forward[ei] << std::endl;
		}
	}

	// Compute log partition function logZ or MAP value
	factor F(0.0);
	for (flist::const_iterator ci = m_roots.begin();
			ci != m_roots.end(); ++ci) {

		factor bel = calc_belief(*ci);
		std::map<size_t, size_t>::iterator mi = m_cluster2var.find(*ci);
		assert(mi != m_cluster2var.end());
		if (m_task == Task::MAR) { // SUM variable
			F += log( bel.sum());
		} else { // MAP variable
			F += log( bel.max() );
		}
	}

	// Partition function or MAP/MMAP value
	m_logz += F.max();

	if (m_debug) {
		std::cout << "Finished forward pass with logZ: " << m_logz << std::endl;
	}

//	size_t C = m_factors.size();
//	for (size_t c = 0; c < C; ++c) {
//		m_factors[c] *= std::exp(m_logz/(double)C);
//	}
}

void ijgp::backward() {

	// Compute and propagate backward messages (bottom-up)
	if (m_debug) {
		std::cout << "Begin backward (bottom-up) pass ..." << std::endl;
	}

	vector<std::pair<findex, findex> >::reverse_iterator ri = m_schedule.rbegin();
	for (; ri != m_schedule.rend(); ++ri ) {

		// Compute backward message m(b->a)
		findex a = (*ri).first;
		findex b = (*ri).second;
		size_t ei = m_edge_indeces[a][b]; // edge index

		// Get the variables to be eliminated
		variable_set VX = m_scopes[b] - m_separators[a][b];

		// compute the belief at b (excluding message a->b)
		factor bel = calc_belief(b, a);
		m_backward[ei] = elim(bel, VX);
		m_backward[ei].normalize();
//		double mx = m_backward[ei].max(); // normalize for stability
//		m_backward[ei] /= mx;

		if (m_debug) {
			std::cout << " - Sending backward msg from " << b << " to " << a << std::endl;
			std::cout << "  - backward msg (" << b << "," << a << "): elim = " << VX << std::endl;
			std::cout << "  -> " << m_backward[ei] << std::endl;
		}
	}

	if (m_debug) {
		std::cout << "Finished backward pass." << std::endl;
	}
}

void ijgp::update() {

	// Compute the marginal (belief) for each variable
	for (vindex v = 0; v < m_gmo.nvar(); ++v) {
		findex c = m_clusters[v][0]; // get a cluster corresp. to current variable
		variable_set vars = m_scopes[c];
		variable VX = m_gmo.var(v);
		//variable_set out = vars - VX;

		factor bel = calc_belief(c);
		m_beliefs[v] = marg(bel, VX);
		if (m_task == Task::MAP) {
			double mx = m_beliefs[v].max();
			m_beliefs[v] /= mx;
		} else {
			m_beliefs[v].normalize();
		}
	}

	// update beliefs (marginals) or compute the MAP/MMAP assignment
	if (m_task == Task::MAR || m_task == Task::PR) {
		for (vindex v = 0; v < m_gmo.nvar(); ++v) {
			findex c = m_clusters[v][0]; // get a cluster corresp. to current variable
			variable_set vars = m_scopes[c];
			variable VX = m_gmo.var(v);

			factor bel = calc_belief(c);
			m_beliefs[v] = marg(bel, VX);
			m_beliefs[v].normalize();
		}
	} else if (m_task == Task::MAP) {
		for (variable_order_t::const_reverse_iterator x = m_order.rbegin();
				x != m_order.rend(); ++x) {

			//variable VX = var(*x);
			findex a = m_clusters[*x][0]; // get source bucket of the variable
			factor bel = incoming(a);

			// condition on previous assignment
			for (variable_order_t::const_reverse_iterator y = m_order.rbegin();
				y != m_order.rend(); ++y) {

				if (*y == *x) break;
				variable VY = var(*y);
				if (m_scopes[a].contains(VY)) {
					bel = bel.condition(VY, m_best_config[*y]);
				}
			}
			m_best_config[*x] = bel.argmax();
		}
	}
}

void ijgp::propagate(size_t nIter, double stopTime, double stopObj) {
	std::cout << "[IJGP] Begin message passing over join graph ..." << std::endl;

	for (size_t iter = 1; iter <= nIter; ++iter) {

		double prevZ = m_logz;
		forward();
		backward();
		update(); // update beliefs

		double dObj = fabs(m_logz - prevZ);
		std::cout << "  logZ: " << std::fixed << std::setw(12)
			<< std::setprecision(MERLIN_PRECISION)
			<< m_logz << " (" << std::scientific
			<< std::setprecision(MERLIN_PRECISION)
			<< std::exp(m_logz) << ") ";
		std::cout << "\td=" << dObj << "\t time="  << std::fixed
			<< std::setprecision(MERLIN_PRECISION)
			<< (timeSystem() - m_start_time)
			<< "\ti=" << iter << std::endl;

		if (dObj < stopObj) {
			break;
		}

		// do at least one iterations
		if (stopTime > 0 && stopTime <= (timeSystem() - m_start_time)) {
			break;
		}
	}
}


} // end namespace



