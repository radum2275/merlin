/*
 * wmb.cpp
 *
 *  Created on: 12 September 2018
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

/// \file wmb.cpp
/// \brief Weighted Mini-Buckets algorithm
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#include "wmb.h"

namespace merlin {

// Write the solution to the output stream
void wmb::write_solution(std::ostream& out, const std::map<size_t, size_t>& evidence,
		const std::map<size_t, size_t>& old2new, const graphical_model& orig,
		const std::set<size_t>& dummies, int output_format) {

	if (output_format == MERLIN_OUTPUT_JSON) {
		out << "{";
		out << " \"algorithm\" : \"wmb\", ";
		out << " \"ibound\" : " << m_ibound << ", ";
		out << " \"iterations\" : " << m_num_iter << ", ";
		switch (m_task) {
		case Task::PR:
			{
				double val = (m_logz + std::log(orig.get_global_const()));
				double prob = std::exp(val);

				out << " \"task\" : \"PR\", ";
				out << " \"value\" : " << std::fixed
					<< std::setprecision(MERLIN_PRECISION)
					<< (m_logz + std::log(orig.get_global_const())) << ", ";
				out << " \"status\" : \"true\" ";

				if (prob == 0.0) { // probability of evidence is 0
					out << " \"status\" : \"false\", ";
					out << " \"message\" : \"Inconsistent evidence or underflow\" ";
				} else {
					out << " \"status\" : \"true\", ";
					out << " \"message\" : \"Consistent evidence\" ";
				}

				break;
			}
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
		case Task::MMAP:
			{
				out << " \"task\" : \"MMAP\", ";
				out << " \"value\" : " << std::fixed
					<< std::setprecision(MERLIN_PRECISION)
					<< (m_logz + std::log(orig.get_global_const())) << ", ";
				out << " \"status\" : \"true\", ";
				out << " \"solution\" : [ ";

				// Evidence variables are a disjoint set from the query variables
				for (vindex i = 0; i < m_query.size(); ++i) {
					vindex j = m_query[i];
					out << "{";
					out << " \"variable\" : " << j << ",";
					assert(m_var_types[j] == true);
					out << " \"value\" : " << m_best_config[j];
					out << "}";
					if (i != m_query.size() - 1) {
						out << ", ";
					}
				}
				out << "] ";
				break;
			}
		}

		out << "}";
	} else if (output_format == MERLIN_OUTPUT_UAI) {
		switch (m_task) {
		case Task::PR:
			{
				double val = m_logz + std::log(orig.get_global_const());
				double prob = std::exp(val);

				out << "PR" << std::endl;
				out << std::fixed << std::setprecision(MERLIN_PRECISION)
					<< (m_logz + std::log(orig.get_global_const()))
					<< " (" << std::scientific << std::setprecision(MERLIN_PRECISION)
					<< std::exp(m_logz + std::log(orig.get_global_const()))
					<< ")" << std::endl;

				if (prob == 0.0) {
					out << "STATUS" << std::endl;
					out << "false: Inconsistent evidence or underflow" << std::endl;
				} else {
					out << "STATUS" << std::endl;
					out << "true: Consistent evidence" << std::endl;
				}

				break;
			}
		case Task::MAR:
			{
				double val = m_logz + std::log(orig.get_global_const());
				double prob = std::exp(val);

				out << "PR" << std::endl;
				out << std::fixed << std::setprecision(MERLIN_PRECISION)
					<< (m_logz + std::log(orig.get_global_const()))
					<< " (" << std::scientific << std::setprecision(MERLIN_PRECISION)
					<< std::exp(m_logz + std::log(orig.get_global_const()))
					<< ")" << std::endl;

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
		case Task::MMAP:
			{
				// evidence variables are a disjoint set from the query variables
				out << "MMAP" << std::endl;
				out << m_query.size();
				for (vindex i = 0; i < m_query.size(); ++i) {
					vindex j = m_query[i];
					assert(m_var_types[j] == true);
					out << " " << m_best_config[j];
				}

				out << std::endl;

				break;
			}
		}
	} else {
		std::string err_msg("Unknown output format.");
		throw std::runtime_error(err_msg);
	}
}

void wmb::init() {

	// Initialize variable types
	m_var_types.resize(m_gmo.nvar(), false); // all SUM
	for (size_t i = 0; i < m_query.size(); ++i) {
		m_var_types[m_query[i]] = true; // set MAP variable
	}

	// Prologue
	std::cout << "[WMB] + i-bound          : " << m_ibound << std::endl;
	std::cout << "[WMB] + iterations       : " << m_num_iter << std::endl;
	std::cout << "[WMB] + inference task   : " << m_task << std::endl;
	if (m_query.empty() == false) {
		std::cout << "+ query vars       : ";
		std::copy(m_query.begin(), m_query.end(), std::ostream_iterator<vindex>(std::cout, " "));
		std::cout << std::endl;
	}
	std::cout << "[WMB] + ordering method  : " << m_order_method << std::endl;
	std::cout << "[WMB] + order iterations : " << m_order_iter << std::endl;
	std::cout << "[WMB] + elimination      : ";

	if (m_order.size() == 0) { // if we need to construct an elimination ordering
		//m_order = m_gmo.order(m_order_method, m_var_types);
		m_order = m_gmo.order2(m_order_method, m_var_types);
//			variable_order_t ord;
//			size_t min_w = 1000000;
//			for (size_t i = 1; i <= m_order_iter; ++i) {
//				ord = m_gmo.order2(m_order_method, m_var_types);
//				size_t w = m_gmo.induced_width(ord);
//				if (w < min_w) {
//					m_order = ord;
//					min_w = w;
//				}
//			}
		m_parents.clear(); // (new elim order => need new pseudotree)
		std::copy(m_order.begin(), m_order.end(),
			std::ostream_iterator<size_t>(std::cout, " "));
	}
	if (m_parents.size() == 0) {     // if we need to construct a pseudo-tree
		m_parents = m_gmo.pseudo_tree(m_order);
	}

	std::cout << std::endl;
	size_t wstar = m_gmo.induced_width(m_order);
	std::cout << "[WMB] + induced width    : " << wstar << std::endl;
	std::cout << "[WMB] + exact inference  : " << (m_ibound >= wstar ? "Yes" : "No") << std::endl;
	std::cout << "[WMB] + ordering time    : " << (timeSystem() - m_start_time) << " seconds" << std::endl;
	if (m_ibound >= wstar) {
		m_num_iter = 1; // exact inference requires 1 iteration over the join-tree
	}

	// Get the factors scopes
	vector<variable_set> fin;
	for (vector<factor>::const_iterator i = m_gmo.get_factors().begin();
			i != m_gmo.get_factors().end(); ++i) {
		fin.push_back((*i).vars());
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

	// First downward pass to initialize the mini-bucket tree and backward messages
	m_clusters.resize(m_order.size());
	for (variable_order_t::const_iterator x = m_order.begin(); x != m_order.end(); ++x) {

		//std::cout << "Eliminating "<<*x << (m_var_types[*x] ? "(MAP)\n" : "(SUM)\n");

		variable VX = var(*x);
		if (*x >= vin.size() || vin[*x].size() == 0)
			continue;  // check that we have some factors over this variable

		flist ids = vin[*x];  // list of factor IDs contained in this bucket

		// Select allocation into mini-buckets
		typedef flist::const_iterator flistIt;
		typedef std::pair<double, Pair> _INS;
		std::multimap<double, Pair> scores;
		std::map<Pair, std::multimap<double, Pair>::iterator> reverseScore;

		// Populate list of pairwise scores for aggregation
		for (flistIt i = ids.begin(); i != ids.end(); ++i) {
			for (flistIt j = ids.begin(); j != i; ++j) {
				double err = score(fin, VX, *i, *j);
				Pair sp(*i, *j);
				reverseScore[sp] = scores.insert(_INS(err, sp)); // save score
			}
			reverseScore[Pair(*i, *i)] = scores.insert(
					_INS(-1, Pair(*i, *i)));       // mark self index at -1
		}

		//   Run through until no more pairs can be aggregated:
		//   Find the best pair (ii,jj) according to the scoring heuristic and join
		//   them as jj; then remove ii and re-score all pairs with jj
		for (;;) {
			std::multimap<double, Pair>::reverse_iterator top =
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
					scores.erase(reverseScore[Pair(ii, *k)]);
				}
				ids /= ii;

				for (flistIt k = ids.begin(); k != ids.end(); ++k) { // updated j; rescore all pairs (j,k)
					if (*k == jj)
						continue;
					double err = score(fin, VX, jj, *k);
					Pair sp(jj, *k);
					scores.erase(reverseScore[sp]);    // change score (i,j)
					reverseScore[sp] = scores.insert(_INS(err, sp));  //
				}
			}
		}

		// Weight for mini-buckets
		double R = (double)ids.size();
		double weight = (m_var_types[*x]) ? infty() : (1.0/R);

		// Eliminate individually each mini-bucket
		vector<findex> alphas;
		int pos = 0;
		for (flistIt i = ids.begin(); i != ids.end(); ++i) {
			//
			// Create new cluster alpha over this set of variables; save function parameters also
			findex alpha = findex(-1);
			alpha = add_factor(factor(fin[*i]));
			alphas.push_back(alpha);
			m_clusters[*x] |= alpha;

			fin[*i] = fin[*i] - VX;

			// add inter clusters edges
			for (flistIt j = New[*i].begin(); j != New[*i].end(); ++j) {
				add_edge(*j, alpha);
				m_schedule.push_back(std::make_pair(*j, alpha));
			}

			// update cluster types
			m_types.push_back(m_var_types[*x]);
			m_weights.push_back(weight);

			// keep track of original factors
			m_originals.push_back(flist());
			m_originals[alpha] |= Orig[*i];
			m_cluster2var[alpha] = *x; // map cluster to variable

			// now incoming nodes to *i is just alpha
			Orig[*i].clear();
			New[*i].clear();
			New[*i] |= alpha;

			// recompute and update adjacency
			insert(vin, *i, fin[*i]);
			++pos;
		}

		//std::cout<<"\n";

	}
	// end for: variable elim order

	// separators and cluster scopes
	size_t C = m_factors.size(), max_clique_size = 0, max_sep_size = 0;
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

	// incoming and outgoing
	m_in.resize(C);
	m_out.resize(C);
	for (vector<std::pair<findex, findex> >::const_iterator i = m_schedule.begin();
			i != m_schedule.end(); ++i) {
		findex from = (*i).first;
		findex to = (*i).second;
		m_in[to] |= from;
		m_out[from] |= to;
	}

	// init the root cluster(s)
	for (size_t i = 0; i < m_out.size(); ++i) {
		if ( m_out[i].empty() )
			m_roots |= i;
	}

	// init forward and backward messages
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

	// init clique potentials
	for (size_t i = 0; i < m_factors.size(); ++i) {
		m_factors[i] = factor(1.0); //get_factor(1.0); // init

		// clique potential
		for (flist::const_iterator j = m_originals[i].begin();
				j != m_originals[i].end(); ++j) {
			m_factors[i] *= m_gmo.get_factor(*j);
		}

	}

	// initialize beliefs (marginals)
	m_logz = 0;
	m_beliefs.clear();
	m_beliefs.resize(m_gmo.nvar(), factor(1.0));
	m_reparam.resize( m_factors.size(), factor(1.0) );
	m_best_config.resize(m_gmo.nvar(), -1);

	// Output the join graph statistics
	std::cout << "[WMB] Created join graph with " << C << " clique factors" << std::endl;
	std::cout << "[WMB] Number of cliques  : " << C << std::endl;
	std::cout << "[WMB] Number of edges    : " << elist.size() << std::endl;
	std::cout << "[WMB] Max clique size    : " << max_clique_size << std::endl;
	std::cout << "[WMB] Max separator size : " << max_sep_size << std::endl;
	std::cout << "[WMB] Finished initialization in " << (timeSystem() - m_start_time) << " seconds" << std::endl;

	if (m_debug) {
		std::cout << "[MERLIN DEBUG]\n";
		std::cout << "[DBG] Join graph with " << m_factors.size() << " clusters and "
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

		// _match list
		std::cout << "[DBG] _MATCH list:" << std::endl;
		for (size_t i = 0; i < m_clusters.size(); ++i) {
			std::cout << "  var " << i << ": ";
			std::copy(m_clusters[i].begin(), m_clusters[i].end(),
					std::ostream_iterator<int>(std::cout, " "));
			std::cout << std::endl;
		}
		std::cout << "[DBG] _WEIGHTS list:" << std::endl;
		for (size_t i = 0; i < m_weights.size(); ++i) {
			std::cout << "  var " << i << ": " << m_weights[i] << std::endl;
		}
		// factors, forward and backward
		std::cout << "[DBG] clique_factors:" << std::endl;
		for (size_t i = 0; i < m_factors.size(); ++i) {
			std::cout << "[" << i << "]: " << m_factors[i] << std::endl;
		}
		std::cout << "[DBG] _forward messages (top-down):" << std::endl;
		for (size_t i = 0; i < m_forward.size(); ++i) {
			std::cout << "(" << i << "): " << m_forward[i] << std::endl;
		}
		std::cout << "[DBG] _backward messages (bottom-up):" << std::endl;
		for (size_t i = 0; i < m_backward.size(); ++i) {
			std::cout << "(" << i << "): " << m_backward[i] << std::endl;
		}
	} // end if debug
}

factor wmb::calc_belief(findex a) {

	factor bel = m_factors[a] * m_reparam[a];

	// forward messages to 'a'
	for (flist::const_iterator ci = m_in[a].begin();
			ci != m_in[a].end(); ++ci) {
		findex p = (*ci);
		size_t j = m_edge_indeces[p][a];
		bel *= m_forward[j];
	}

	// backward message to 'a'
	for (flist::const_iterator ci = m_out[a].begin();
			ci != m_out[a].end(); ++ci) {
		findex p = (*ci);
		size_t j = m_edge_indeces[a][p];
		bel *= m_backward[j];
	}

	return bel;
}

factor wmb::incoming(findex a, size_t i) {

	factor bel = m_factors[a] * m_reparam[a];

	// forward messages to 'a'
	for (flist::const_iterator ci = m_in[a].begin();
			ci != m_in[a].end(); ++ci) {
		findex p = (*ci);
		size_t j = m_edge_indeces[p][a];
		bel *= m_forward[j];
		//_norm[i] += _norm[j];
	}

	return bel;
}

factor wmb::incoming(findex a) {

	factor bel = m_factors[a] * m_reparam[a];

	// forward messages to 'a'
	for (flist::const_iterator ci = m_in[a].begin();
			ci != m_in[a].end(); ++ci) {
		findex p = (*ci);
		size_t j = m_edge_indeces[p][a];
		bel *= m_forward[j];
	}

	return bel;
}

void wmb::forward(double step) {

	if (m_debug) {
		std::cout << "Begin forward (top-down) pass ..." << std::endl;
	}

	m_logz = 0;
	for (variable_order_t::const_iterator x = m_order.begin(); x != m_order.end(); ++x) {

		if (m_debug) {
			std::cout << " - Eliminating " << *x
				<< (m_var_types[*x] ? " (MAP)\n" : " (SUM)\n");
		}

		// Moment-match the clusters of this bucket
		match_clusters(*x, step);

		// Generate forward messages from each of the clusters corresp. to x
		variable VX = var(*x);
		for (flist::const_iterator it = m_clusters[*x].begin();
				it != m_clusters[*x].end(); ++it) {
			findex a = (*it);
			if ( m_out[a].size() > 0 ) {
				findex b = *(m_out[a].begin());
				size_t ei = m_edge_indeces[a][b];

				factor tmp = incoming(a, ei);
				if (m_var_types[*x] == false) { // SUM
					m_forward[ei] = tmp.sum_power(VX, 1.0/m_weights[a]);
				} else { // MAX
					m_forward[ei] = tmp.max(VX);
				}

				// normalize for numerical stability
				double mx = m_forward[ei].max();
				m_forward[ei] /= mx;
				m_logz += std::log(mx);

				if (m_debug) {
					std::cout << "  forward msg (" << a << "," << b << "): elim = " << VX << " -> ";
					std::cout << m_forward[ei] << std::endl;
				}
			} // end if
		} // end for
	} // end for

	// Compute log partition function logZ or MAP/MMAP value
	factor F(0.0);
	for (flist::const_iterator ci = m_roots.begin();
			ci != m_roots.end(); ++ci) {

		factor bel = calc_belief(*ci);
		std::map<size_t, size_t>::iterator mi = m_cluster2var.find(*ci);
		assert(mi != m_cluster2var.end());
		size_t v = mi->second;
		if (m_var_types[v] == false) { // SUM variable
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
}

void wmb::backward(size_t iter) {

	if (m_debug) {
		std::cout << "Begin backward (bottom-up) pass ..." << std::endl;
	}

	// update backward messages
	vector<std::pair<findex, findex> >::reverse_iterator ri = m_schedule.rbegin();
	for (; ri != m_schedule.rend(); ++ri ) {

		// compute backward message m(b->a)
		findex a = (*ri).first;
		findex b = (*ri).second;
		size_t i = m_edge_indeces[a][b]; // edge index

		variable_set VX = m_scopes[b] - m_separators[a][b];

		if (m_debug) {
			std::cout << " - Sending backward msg from " << a << " to " << b << std::endl;
		}

		// compute the belief at b
		factor bel = calc_belief(b);

		if (m_types[b] == false && m_types[a] == false) { // SUM-SUM

			bel ^= 1.0/m_weights[b];
			bel /= (m_forward[i]^(1.0/m_weights[a])); // divide out m(a->b)

			//_backward[i] = elim(bel, VX, 1);
			m_backward[i] = bel.sum(VX);
			m_backward[i] ^= (m_weights[a]);

		} else if (m_types[b] == true && m_types[a] == true) { // MAX-MAX

			bel /= m_forward[i]; // divide out m(a->b)
			m_backward[i] = bel.max(VX);

		} else if (m_types[b] == true && m_types[a] == false) { // MAX-SUM

			bel = bel.sigma(iter); // the sigma operator that focuses on max
			bel /= (m_forward[i]^(1.0/m_weights[a])); // divide out m(a->b)

			m_backward[i] = bel.sum(VX);
			m_backward[i] ^= (m_weights[a]);

		} else {
			assert(false); // cannot reach this case!!
		}

		// normalize for numerical stability
		double mx = m_backward[i].max();
		m_backward[i] /= mx;

		if (m_debug) {
			std::cout << "  backward msg (" << b << "," << a << "): elim = " << VX << " -> ";
			std::cout << m_backward[i] << std::endl;
		}

	}

	if (m_debug) {
		std::cout << "Finished backward (bottom-up) pass." << std::endl;
	}
}

void wmb::match_clusters(size_t x, double step) {

	if (m_clusters[x].size() <= 1) {
		return; // no matching
	}

	variable VX = var(x);
	if (m_var_types[x] == true) { // max marginals matching

		size_t R = m_clusters[x].size();
		vector<factor> ftmp(R); // compute geometric mean

		variable_set var;
		var |= VX; // on mutual variable (bucket variable)
		factor fmatch(var,1.0);

		size_t i = 0;
		for (flist::const_iterator it = m_clusters[x].begin();
				it != m_clusters[x].end(); ++it, i++) {

			findex a = (*it);
			factor bel = calc_belief(a);
			ftmp[i] = bel.maxmarginal(var); // max-marginal
			fmatch *= ftmp[i];
		}

		i = 0;
		fmatch ^= (1.0/R); // and match each bucket to it
		for (flist::const_iterator it = m_clusters[x].begin();
				it != m_clusters[x].end(); ++it, ++i) {
			findex a = (*it);
			m_reparam[a] *= (fmatch/ftmp[i]);
		}

	} else { // weighted marginals matching

		size_t R = m_clusters[x].size();
		vector<factor> ftmp(R);   // compute geometric mean

		variable_set var;
		var |= VX; // on mutual variable (bucket variable)
		factor fmatch(var,1.0);
		size_t i = 0;

		for (flist::const_iterator it = m_clusters[x].begin();
				it != m_clusters[x].end(); ++it, i++) {

			findex a = (*it);
			factor bel = calc_belief(a);
			bel ^= (1.0/m_weights[a]);
			ftmp[i] = bel.marginal(var);
			fmatch *= (ftmp[i] ^ m_weights[a]);
		}

		i = 0;
		for (flist::const_iterator it = m_clusters[x].begin();
				it != m_clusters[x].end(); ++it, ++i) {
			findex a = (*it);
			m_reparam[a] *= ((fmatch/ftmp[i])^(step*m_weights[a]));
		}
	}
}

void wmb::tighten(size_t nIter, double stopTime, double stopObj) {
	std::cout << "[WMB] Begin message passing over join graph ..." << std::endl;

	double minZ = infty();
	for (size_t iter = 1; iter <= nIter; ++iter) {
		double step = 1.0/(double)iter;
		double prevZ = m_logz;

		forward(step);
		backward(iter);
		update();

		// keep track of tightest upper bound
		if (m_logz < minZ) {
			minZ = m_logz;
		}

		double dObj = fabs(m_logz - prevZ);
		std::cout << "  logZ: " << std::fixed
				<< std::setw(12) << std::setprecision(MERLIN_PRECISION)
				<< m_logz << " (" << std::scientific
				<< std::setprecision(MERLIN_PRECISION)
				<< std::exp(m_logz) << ") ";
		std::cout << "\td=" << dObj << "\t time="  << std::fixed
			<< std::setprecision(MERLIN_PRECISION) << (timeSystem() - m_start_time)
			<< "\ti=" << iter << std::endl;

		if (dObj < stopObj) {
			break;
		}

		// do at least one iterations
		if (stopTime > 0 && stopTime <= (timeSystem() - m_start_time)) {
			break;
		}
	} // end for

	m_logz = minZ; // keep tightest upper bound
}

void wmb::update() {

	// update beliefs (marginals) or compute the MAP/MMAP assignment
	if (m_task == Task::MAR || m_task == Task::PR) {
		for (vindex v = 0; v < m_gmo.nvar(); ++v) {
			findex c = m_clusters[v][0]; // get a cluster corresp. to current variable
			double w = m_weights[c];
			variable_set vars = m_scopes[c];
			variable VX = m_gmo.var(v);
			variable_set out = vars - VX;

			factor bel = calc_belief(c);
			m_beliefs[v] = marg(bel, VX, w);
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
	} else if (m_task == Task::MMAP) {
		for (variable_order_t::const_reverse_iterator x = m_order.rbegin();
				x != m_order.rend(); ++x) {

			if (m_var_types[*x] == false) break; // stop at first SUM variable
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

void wmb::run() {

	// Start the timer and store it
	m_start_time = timeSystem();

	init();
	tighten(m_num_iter);

	// Output solution (UAI output format)
	std::cout << "[WMB] Converged after " << m_num_iter << " iterations in "
			<< (timeSystem() - m_start_time) << " seconds" << std::endl;

	switch (m_task) {
	case Task::PR:
	case Task::MAR:
		{
			std::cout << "PR" << std::endl;
			std::cout << std::fixed << std::setprecision(MERLIN_PRECISION)
				<< m_logz << " (" << std::scientific
				<< std::setprecision(MERLIN_PRECISION)
				<< std::exp(m_logz)	<< ")" << std::endl;

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
			std::cout << "Final Upper Bound is " << std::setprecision(MERLIN_PRECISION)
				<< m_logz << " (" << std::scientific << std::setprecision(MERLIN_PRECISION)
				<< std::exp(m_logz) << ")" << std::endl;
			std::cout << "Final Lower Bound is " << std::setprecision(MERLIN_PRECISION)
				<< m_lb << " (" << std::scientific << std::setprecision(MERLIN_PRECISION)
				<< std::exp(m_lb) << ")" << std::endl;
			std::cout << "MAP" << std::endl;
			std::cout << m_gmo.nvar();
			for (vindex v = 0; v < m_gmo.nvar(); ++v) {
				std::cout << " " << m_best_config[v];
			}
			std::cout << std::endl;

			break;
		}
	case Task::MMAP:
		{
			std::cout << "Final Upper Bound is " << std::setprecision(MERLIN_PRECISION)
				<< m_logz << " (" << std::scientific << std::setprecision(MERLIN_PRECISION)
				<< std::exp(m_logz) << ")" << std::endl;
			std::cout << "MMAP" << std::endl;
			std::cout << m_query.size();
			for (vindex v = 0; v < m_gmo.nvar(); ++v) {
				if (m_var_types[v] == true)
					std::cout << " " << m_best_config[v];
			}
			std::cout << std::endl;
			break;
		}
	} // end switch

}


} // end namespace



