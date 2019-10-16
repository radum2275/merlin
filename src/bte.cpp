/*
 * bte.cpp
 *
 *  Created on: 12 September 2018
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

/// \file bte.cpp
/// \brief Bucket Tree Elimination algorithm
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#include "bte.h"

namespace merlin {

///
/// \brief Write the solution to the output stream
///
void bte::write_solution(std::ostream& out, const std::map<size_t, size_t>& evidence,
		const std::map<size_t, size_t>& old2new, const graphical_model& orig,
		const std::set<size_t>& dummies, int output_format) {

	if (output_format == MERLIN_OUTPUT_JSON) {
		out << "{";
		out << " \"algorithm\" : \"bte\", ";
		switch (m_task) {
		case Task::PR:
			{
				double val = m_logz + std::log(orig.get_global_const());
				double prob = std::exp(val);

				out << " \"task\" : \"PR\", ";
				out << " \"value\" : " << std::fixed
					<< std::setprecision(MERLIN_PRECISION)
					<< val << ", ";
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
						variable v = orig.var(i);
						if (dummies.find(i) != dummies.end()) {
							continue; // skip dummy variables from output
						}

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
		case Task::MAR:
			{
				double val = m_logz + std::log(orig.get_global_const());
				double prob = std::exp(val);

				out << "PR" << std::endl;
				out << std::fixed << std::setprecision(MERLIN_PRECISION)
					<< (m_logz + std::log(orig.get_global_const())) << " ("
					<< std::scientific << std::setprecision(MERLIN_PRECISION)
					<< std::exp(m_logz + std::log(orig.get_global_const()))
					<< ")" << std::endl;

				if (prob == 0.0) { // probability of evidence is 0
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

///
/// \brief Run the bucket-elimination algorithm
///
void bte::run() {

	// Start the timer and store it
	m_start_time = timeSystem();

	init();
	propagate();

	// Output solution (UAI output format)
	std::cout << "[BTE] Finished in " << (timeSystem() - m_start_time)
			<< " seconds" << std::endl;


	switch (m_task) {
	case Task::PR:	// partition function
		{
			std::cout << "PR" << std::endl;
			std::cout << std::fixed << std::setprecision(MERLIN_PRECISION)
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

			break;
		}
	case Task::MAR:	// marginals
		{
			std::cout << "PR" << std::endl;
			std::cout << std::fixed << std::setprecision(MERLIN_PRECISION)
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
				for (size_t k = 0; k < VX.states(); ++k) {
					std::cout << " " << std::fixed
						<< std::setprecision(MERLIN_PRECISION)
						<< belief(VX)[k];
				}
			}
			std::cout << std::endl;

			break;
		}
	case Task::MAP: // MAP inference
		{
			std::cout << "Value" << std::endl;
			std::cout << std::fixed << std::setprecision(MERLIN_PRECISION)
				<< m_logz << " (" << std::scientific
				<< std::setprecision(MERLIN_PRECISION)
				<< std::exp(m_logz) << ")" << std::endl;

			std::cout << "MAP" << std::endl;
			std::cout << m_gmo.nvar();
			for (vindex v = 0; v < m_gmo.nvar(); ++v) {
				std::cout << " " << m_best_config[v];
			}
			std::cout << std::endl;

			break;
		}
	case Task::MMAP: // Marginal MAP inference
		{
			std::cout << "Value" << std::endl;
			std::cout << std::fixed << std::setprecision(MERLIN_PRECISION)
				<< m_logz << " (" << std::scientific
				<< std::setprecision(MERLIN_PRECISION)
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
	}

}

///
/// \brief Initialize the bucket-elimination algorithm
///
void bte::init() {

	// Initialize variable types
	m_var_types.resize(m_gmo.nvar(), false); // all SUM
	for (size_t i = 0; i < m_query.size(); ++i) {
		vindex qvar = m_query[i];
		m_var_types[qvar] = true; // set MAP variable
	}

	// Prologue
	std::cout << "[BTE] + inference task   : " << m_task << std::endl;
	if (m_query.empty() == false) {
		std::cout << "+ query vars       : ";
		std::copy(m_query.begin(), m_query.end(), std::ostream_iterator<vindex>(std::cout, " "));
		std::cout << std::endl;
	}
	std::cout << "[BTE] + ordering method  : " << m_order_method << std::endl;
	std::cout << "[BTE] + elimination      : ";

	if (m_order.size() == 0) { // if we need to construct an elimination ordering
		m_order = m_gmo.order(m_order_method, m_var_types);
		m_parents.clear(); // (new elim order => need new pseudotree)
		std::copy(m_order.begin(), m_order.end(),
				std::ostream_iterator<size_t>(std::cout, " "));
	}
	if (m_parents.size() == 0) { // if we need to construct a pseudo tree
		m_parents = m_gmo.pseudo_tree(m_order);
	}

	// Calculate the induced width of the elimination ordering
	std::cout << std::endl;
	size_t wstar = m_gmo.induced_width(m_order);
	std::cout << "[BTE] + induced width    : " << wstar << std::endl;
	std::cout << "[BTE] + exact inference  : Yes" << std::endl;
	std::cout << "[BTE] + ordering time    : " << (timeSystem() - m_start_time) << " seconds" << std::endl;

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

	vector<flist> Orig(m_gmo.num_factors()); // origination info: which original factors are
	for (size_t i = 0; i < Orig.size(); ++i)
		Orig[i] |= i;    					// included for the first time, and which newly
	vector<flist> New(m_gmo.num_factors()); // created clusters feed into this cluster

	// First downward pass to initialize the bucket tree and backward messages
	if (m_debug) {
		std::cout << "[BTE] Initializing bucket-tree ... " << std::endl;
	}

	m_clusters.resize(m_order.size());
	for (variable_order_t::const_iterator x = m_order.begin(); x != m_order.end(); ++x) {

		if (m_debug) {
			std::cout << "  - create bucket/cluster for var " << *x
					<< (m_var_types[*x] ? " (MAP)\n" : " (SUM)\n");
		}

		// Get the current variable to process
		variable VX = var(*x);
		if (*x >= vin.size() || vin[*x].size() == 0)
			continue;  // check that we have some factors over this variable

		flist ids = vin[*x];  // list of factor IDs contained in this bucket

		if (m_debug) {
			std::cout << "  - factors in this bucket: " << ids.size() << std::endl;
			for (flist::const_iterator i = ids.begin(); i != ids.end(); ++i) {
				std::cout << "   original factor id " << *i << " : " << fin[*i] << " --> " << m_gmo.get_factor(*i) << std::endl;
			}
		}

		// Merge all factor scopes in this bucket into a single scope
		vector<size_t> temp;
		typedef flist::const_iterator flistIt;
		size_t jj = *(ids.begin());
		for (flistIt i = ids.begin(); i != ids.end(); ++i) {
			size_t ii = *i;
			if (ii == jj) continue;
			fin[jj] |= fin[ii];     // combine into j
			erase(vin, ii, fin[ii]);
			fin[ii] = variable_set();  //   & remove i

			Orig[jj] |= Orig[ii];
			Orig[ii].clear(); // keep track of list of original factors in this cluster
			New[jj] |= New[ii];
			New[ii].clear(); //  list of new "message" clusters incoming to this cluster

			temp.push_back(ii);
		}

		for (size_t i = 0; i < temp.size(); ++i) {
			size_t ii = temp[i];
			ids /= ii;
		}

		// Sanity checks
		assert(ids.size() == 1);
		if (m_debug) {
			std::cout << "  After merging: " << ids.size() << std::endl;
			for (flist::const_iterator i = ids.begin(); i != ids.end(); ++i) {
				std::cout << "  Factor id " << *i << std::endl;
				std::cout << "  Scope: " << fin[*i] << std::endl;
			}
		}

		// Eliminate the bucket variable
		vector<findex> alphas;
		for (flistIt i = ids.begin(); i != ids.end(); ++i) {
			//
			// Create new cluster alpha over this set of variables; save function parameters also
			findex alpha = findex(-1);
			alpha = add_factor(factor(fin[*i])); // add the clique factor
			alphas.push_back(alpha);

			fin[*i] = fin[*i] - VX;

			// Add inter clusters edges
			for (flistIt j = New[*i].begin(); j != New[*i].end(); ++j) {
				add_edge(*j, alpha);
				m_schedule.push_back(std::make_pair(*j, alpha));
			}

			// Keep track of original factors
			m_originals.push_back(flist());
			m_originals[alpha] |= Orig[*i];
			m_clusters[*x] |= alpha; // map variable to bucket/cluster
			m_cluster2var[alpha] = *x; // map cluster to variable

			// Now incoming nodes to *i is just alpha
			Orig[*i].clear();
			New[*i].clear();
			New[*i] |= alpha;

			// Recompute and update adjacency
			insert(vin, *i, fin[*i]);
		}
	}
	// end for: variable elim order
	if (m_debug) {
		std::cout << "  - number of clique factors is: " << m_factors.size() << std::endl;
		std::cout << "[BTE] Done initializing the bucket-tree." << std::endl;
	}

	std::cout << "[BTE] Created bucket-tree with " << m_factors.size()
			<< " clique factors" << std::endl;

	// Set the separators and cluster scopes
	size_t max_sep_size = 0, max_clique_size = 0;
	size_t C = m_factors.size(); // number of clique factors
	m_separators.resize(C);
	for (size_t i = 0; i < C; ++i) m_separators[i].resize(C);
	m_scopes.resize(C); // clique scopes
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

	// Set the incoming and outgoing lists for each cluster/clique
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

	// Init forward and backward messages
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

	// Init clique potentials (un-normalized)
	for (size_t c = 0; c < m_factors.size(); ++c) {
		m_factors[c] = factor(1.0); // init

		// clique potential
		for (flist::const_iterator j = m_originals[c].begin();
				j != m_originals[c].end(); ++j) {
			m_factors[c] *= m_gmo.get_factor(*j);
		}
	}

	// Initialize beliefs (marginals)
	m_logz = 0;
	m_beliefs.clear();
	m_beliefs.resize(m_gmo.nvar(), factor(1.0));
	m_best_config.resize(m_gmo.nvar(), -1);

	// Output the bucket tree statistics
	double elapsed = (timeSystem() - m_start_time);
	std::cout << "[BTE] Number of cliques  : " << C << std::endl;
	std::cout << "[BTE] Number of edges    : " << elist.size() << std::endl;
	std::cout << "[BTE] Max clique size    : " << max_clique_size << std::endl;
	std::cout << "[BTE] Max separator size : " << max_sep_size << std::endl;
	std::cout << "[BTE] Finished initialization in " << elapsed << " seconds" << std::endl;

	if (m_debug) {
		std::cout << std::endl;
		std::cout << "[MERLIN DEBUG]\n";
		std::cout << "[DBG] Bucket-tree with " << m_factors.size() << " clusters and "
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
	} // end if debug
}

///
/// \brief Forward (top-down) message passing along the edge of the bucket tree.
///
void bte::forward() {

	if (m_debug) {
		std::cout << "Begin forward (top-down) pass ..." << std::endl;
	}

	m_logz = 0; // reset the log parition function
	double timestamp = timeSystem();
	for (variable_order_t::const_iterator x = m_order.begin(); x != m_order.end(); ++x) {

		if (m_debug) {
			std::cout << " - Eliminating var " << *x
					<< (m_var_types[*x] ? " (MAP)\n" : " (SUM)\n");
		}

		// Generate forward messages
		variable VX = var(*x);
		findex a = m_clusters[*x][0]; // get source bucket of the variable
		if (m_out[a].size() > 0) {
			findex b = *(m_out[a].begin()); // destination bucket
			size_t ei = m_edge_indeces[a][b]; // edge index (for message)
			if (m_var_types[*x] == false) { // SUM variable
				m_forward[ei] = incoming(a, ei).sum(VX);
			} else { // MAP variable
				m_forward[ei] = incoming(a, ei).max(VX);
			}

			if (m_debug) {
				std::cout << "  forward msg (" << a << "," << b << "): elim = " << VX << " -> ";
				std::cout << m_forward[ei] << std::endl;
			}
		}
	} // done

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

	std::cout << "[BTE] Finished forward pass in " << (timeSystem() - timestamp)
			<< " seconds" << std::endl;
}

///
/// \brief Backward (bottom-up) message passing along the edges of the bucket tree.
///
void bte::backward() {

	// Calculate 'backward' messages
	if (m_debug) {
		std::cout << "Begin backward (bottom-up) pass ..." << std::endl;
	}

	vector<std::pair<findex, findex> >::reverse_iterator ri = m_schedule.rbegin();
	double timestamp = timeSystem();
	for (; ri != m_schedule.rend(); ++ri ) {

		// Select backward message m(b->a)
		findex a = (*ri).first;
		findex b = (*ri).second;
		size_t i = m_edge_indeces[a][b]; // edge index

		// Select variables to eliminate
		variable_set VX = m_scopes[b] - m_separators[a][b];

		if (m_debug) {
			std::cout << " - Sending backward msg from "
				<< a << " to " << b << std::endl;
		}

		// Compute the belief at b
		factor bel = calc_belief(b);
		bel /= m_forward[i]; // discard the forward message on the edge
		if (m_task == Task::PR || m_task == Task::MAR) { // marginals
			m_backward[i] = bel.sum(VX);
		} else if (m_task == Task::MAP) { // MAP inference
			m_backward[i] = bel.max(VX);
		} else { // Marginal MAP inference
			for (variable_set::const_iterator vi = VX.begin();
					vi != VX.end(); ++vi) {
				size_t v = vi->label();
				if (m_var_types[v] == false) { // SUM variable
					bel = bel.sum(*vi);
				} else { // MAP variable
					bel = bel.max(*vi);
				}
			}
			m_backward[i] = bel;
		}

		if (m_debug) {
			std::cout << "  - backward msg (" << b << "," << a << "): elim = " << VX << std::endl;
			std::cout << "  -> " << m_backward[i] << std::endl;
		}
	}

	std::cout << "[BTE] Finished backward pass in " << (timeSystem() - timestamp)
			<< " seconds" << std::endl;
}

///
/// \brief Propagate the messages along the edges of the bucket tree.
///
void bte::propagate() {

	// Top-down 'forward' message passing
	forward();

	if (m_task == Task::PR) {
		return; // stop here to compute the log partition function
	}

	// Bottom-up 'backward' message passing
	backward();

	// Update beliefs/marginals
	update();
}

///
/// \brief Update the beliefs (marginals or max-marginals)
///
void bte::update() {

	// update beliefs (marginals) or compute the MAP/MMAP assignment
	if (m_task == Task::MAR) {
		for (vindex v = 0; v < m_gmo.nvar(); ++v) {
			findex c = m_clusters[v][0]; // get a cluster corresp. to current variable
			variable_set vars = m_scopes[c];
			variable VX = m_gmo.var(v);
			variable_set out = vars - VX;

			factor bel = calc_belief(c);
			m_beliefs[v] = marg(bel, VX);
			m_beliefs[v].normalize();
			//m_beliefs[v] /= std::exp(m_logz); // normalize by logZ
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

factor bte::calc_belief(findex a) {

	factor bel = m_factors[a]; // original factors in the cluster

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

factor bte::incoming(findex a, size_t i) {

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

factor bte::incoming(findex a) {

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

} // end namespace


