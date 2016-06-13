/*
 * be.h
 *
 *  Created on: 25 May 2016
 *      Author: radu
 */

/// \file be.h
/// \brief Bucket (Tree) Elimination algorithm
/// \author Radu Marinescu

#ifndef IBM_MERLIN_BE_H_
#define IBM_MERLIN_BE_H_

#include "graphical_model.h"

namespace merlin {

/**
 * Bucket (Tree) Elimination (BE)
 *
 * Tasks supported: PR, MAR, MAP, MMAP
 *
 * Bucket Tree Elimination is an exact inference algorithm that supports the
 * main inference tasks in graphical models (PR, MAR, MAP, MMAP). It belongs to
 * the class on join-tree elimination algorithms and is driven by a bucket tree
 * structure. Its complexity is time and space exponential in the induced width.
 * Note that this is meant to be a reference implementation and therefore is
 * applicable to relatively small models having relatively small induced widths.
 */
class be: public graphical_model, public algorithm {
public:
	typedef graphical_model::findex findex;        ///< Factor index
	typedef graphical_model::vindex vindex;        ///< Variable index
	typedef graphical_model::flist flist;          ///< Collection of factor indices

public:

	///
	/// \brief Default constructor.
	///
	be() : graphical_model() {
		set_properties();
		m_debug = false;
	}

	///
	/// \brief Constructor with a graphical model.
	///
	be(const graphical_model& gm) : graphical_model(gm), m_gmo(gm) {
		clear_factors();
		set_properties();
		m_debug = false;
	}

	///
	/// \brief Clone the algorithm.
	/// \return the pointer to the new object containing the cloned algorithm.
	///
	virtual be* clone() const {
		be* gm = new be(*this);
		return gm;
	}

	///
	double ub() const {
		return m_log_z;
	}
	double lb() const {
		throw std::runtime_error("Not implemented");
	}
	std::vector<size_t> best_config() const {
		return m_best_config;
	}

	double logZ() const {
		return m_log_z;
	}
	double logZub() const {
		return m_log_z;
	}
	double logZlb() const {
		return m_log_z;
	}

	const factor& belief(size_t f) const {
		return m_beliefs[f];
	}
	const factor& belief(variable v) const {
		return m_beliefs[v];
	}
	const factor& belief(variable_set vs) const {
		throw std::runtime_error("Not implemented");
	}
	const vector<factor>& beliefs() const {
		return m_beliefs;
	}

	///
	/// \brief Get the original graphical model.
	///
	const graphical_model& get_gm_orig() const {
		return m_gmo;
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

		switch (m_task) {
		case Task::PR:
			{
				out << "PR" << std::endl;
				out << std::fixed << std::setprecision(6)
					<< (m_log_z + std::log(orig.get_global_const())) << " (" << std::scientific << std::setprecision(6)
					<< std::exp(m_log_z + std::log(orig.get_global_const())) << ")" << std::endl;

				break;
			}
		case Task::MAR:
			{
				out << "PR" << std::endl;
				out << std::fixed << std::setprecision(6)
					<< (m_log_z + std::log(orig.get_global_const())) << " (" << std::scientific << std::setprecision(6)
					<< std::exp(m_log_z + std::log(orig.get_global_const())) << ")" << std::endl;

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
	}

	///
	/// \brief Run the bucket (tree) elimination algorithm.
	///
	virtual void run() {
		m_debug = true;
		init();
		propagate();

		// Output solution (UAI output format)
		std::cout << "Finished in " << (timeSystem() - m_start_time)
				<< " seconds" << std::endl;


		switch (m_task) {
		case Task::PR:	// partition function
			{
				std::cout << "PR" << std::endl;
				std::cout << std::fixed << std::setprecision(6)
					<< m_log_z << " (" << std::scientific << std::setprecision(6)
					<< std::exp(m_log_z) << ")" << std::endl;
				break;
			}
		case Task::MAR:	// marginals
			{
				std::cout << "PR" << std::endl;
				std::cout << std::fixed << std::setprecision(6)
					<< m_log_z << " (" << std::scientific << std::setprecision(6)
					<< std::exp(m_log_z) << ")" << std::endl;
				std::cout << "MAR" << std::endl;
				std::cout << m_gmo.nvar();
				for (vindex v = 0; v < m_gmo.nvar(); ++v) {
					variable VX = m_gmo.var(v);
					std::cout << " " << VX.states();
					for (size_t j = 0; j < VX.states(); ++j)
						std::cout << " " << std::fixed << std::setprecision(6) << belief(VX)[j];
				}
				std::cout << std::endl;

				break;
			}
		case Task::MAP: // MAP inference
			{
				std::cout << "Value" << std::endl;
				std::cout << std::fixed << std::setprecision(6)
					<< m_log_z << " (" << std::scientific << std::setprecision(6)
					<< std::exp(m_log_z) << ")" << std::endl;

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
				std::cout << std::fixed << std::setprecision(6)
					<< m_log_z << " (" << std::scientific << std::setprecision(6)
					<< std::exp(m_log_z) << ")" << std::endl;

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
	/// \brief Inference tasks supported.
	///
	MER_ENUM( Task, PR,MAR,MAP,MMAP );

	///
	/// \brief Properties of the algorithm
	///
	MER_ENUM( Property , Order,Task );


	// Setting properties (directly or through property string):

	///
	/// \brief Set the variable types.
	///
	void set_var_types(const vector<bool>& var_types) {
		m_var_types = var_types;
	}

	///
	/// \brief Get the variable types.
	///
	const vector<bool>& get_var_types() const {
		return m_var_types;
	}

	///
	/// \brief Set the variable order.
	///
	void set_order(const variable_order_t& ord) {
		m_order = ord;
	}

	///
	/// \brief Set the variable order method.
	///
	void set_order_method(graphical_model::OrderMethod method) {
		m_order.clear();
		m_order_method = method;
	}

	///
	/// \brief Get the variable order.
	///
	const variable_order_t& get_order() {
		return m_order;
	}

	///
	/// \brief Get the pseudo tree.
	///
	const std::vector<vindex>& get_pseudo_tree() {
		return m_parents;
	}

	///
	/// \brief Set the pseudo tree.
	///
	void set_pseudo_tree(const vector<vindex>& p) {
		m_parents = p;
	}

	///
	/// \brief Set the query (MAP) variables.
	///
	void set_query(const std::vector<vindex>& q) {
		m_query = q;
	}

	///
	/// \brief Get the query (MAP) variables.
	///
	const std::vector<vindex>& get_query() {
		return m_query;
	}

	///
	/// \brief Set the graphical model.
	///
	void set_graphical_model(const graphical_model& gm) {
		m_gmo = gm;
	}

	///
	/// \brief Set the graphical model from a list of factors.
	///
	void set_graphical_model(const vector<factor>& fs) {
		m_gmo = graphical_model(fs);
	}

	///
	/// \brief Set the properties of the algorithm.
	/// \param opt 	The string containing comma separated property value pairs
	///
	virtual void set_properties(std::string opt = std::string()) {
		if (opt.length() == 0) {
			set_properties("Order=MinFill,Task=PR");
			return;
		}
		std::vector<std::string> strs = merlin::split(opt, ',');
		for (size_t i = 0; i < strs.size(); ++i) {
			std::vector<std::string> asgn = merlin::split(strs[i], '=');
			switch (Property(asgn[0].c_str())) {
			case Property::Order:
				m_order.clear();
				m_parents.clear();
				m_order_method = graphical_model::OrderMethod(asgn[1].c_str());
				break;
			case Property::Task:
				m_task = Task(asgn[1].c_str());
				break;
			default:
				break;
			}
		}
	}

	///
	/// \brief Eliminate a set of variables from a factor.
	/// \param F 	The reference of the factor to eliminate from
	///	\param vs 	The set of variables to be eliminated
	/// \param w 	The weight of the weighted elimination operator
	/// \return the factor resulted from eliminating the set of variables.
	///
	factor elim(const factor& F, const variable_set& vs, const double w) {
		return F.sum_power(vs, w);
	}

	///
	/// \brief Compute the weighted marginal over a set of variables.
	/// \param F 	The reference of the factor to marginalize over
	/// \param vs 	The set of variables representing the scope of the marginal
	/// \param w 	The weight of the weighted elimination operator
	/// \return the factor representing the weighted marginal over the set of variables.
	///
	factor marg(const factor& F, const variable_set& vs) {
		return F.marginal(vs);
	}

	///
	/// \brief Initialize the bucket (tree) elimination algorithm.
	///
	void init() {

		// Start the timer and store it
		m_start_time = timeSystem();

		// Initialize variable types
		m_var_types.resize(m_gmo.nvar(), false); // all SUM
		for (size_t i = 0; i < m_query.size(); ++i) {
			vindex qvar = m_query[i];
			m_var_types[qvar] = true; // set MAP variable
		}

		// Prologue
		std::cout << VERSIONINFO << std::endl << COPYRIGHT << std::endl;
		std::cout << "Initialize inference engine ..." << std::endl;
		std::cout << "+ tasks supported  : PR, MAR, MAP, MMAP" << std::endl;
		std::cout << "+ algorithm        : " << "BE" << std::endl;
		std::cout << "+ inference task   : " << m_task << std::endl;
		if (m_query.empty() == false) {
			std::cout << "+ query vars       : ";
			std::copy(m_query.begin(), m_query.end(), std::ostream_iterator<vindex>(std::cout, " "));
			std::cout << std::endl;
		}
		std::cout << "+ ordering heur    : " << m_order_method << std::endl;
		std::cout << "+ elimination      : ";

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
		std::cout << "+ induced width    : " << wstar << std::endl;
		std::cout << "+ exact inference  : Yes" << std::endl;

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
		std::cout << "Initializing bucket-tree ... " << std::endl;
		std::cout << "  - initial number of clique factors is: " << m_factors.size() << std::endl;
		m_clusters.resize(m_order.size());
		for (variable_order_t::const_iterator x = m_order.begin(); x != m_order.end(); ++x) {

			std::cout << "  - create bucket/cluster for var " << *x
					<< (m_var_types[*x] ? " (MAP)\n" : " (SUM)\n");

			// Get the current variable to process
			variable VX = var(*x);
			if (*x >= vin.size() || vin[*x].size() == 0)
				continue;  // check that we have some factors over this variable

			flist ids = vin[*x];  // list of factor IDs contained in this bucket

			if (m_debug) {
				std::cout << "  - factors in this bucket: " << ids.size() << std::endl;
				for (flist::const_iterator i = ids.begin(); i != ids.end(); ++i) {
					std::cout << "   factor id " << *i << " : " << fin[*i] << std::endl;
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
		std::cout << "  - final number of clique factors is: " << m_factors.size() << std::endl;
		std::cout << "Done initializing the bucket-tree." << std::endl;

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
		m_log_z = 0;
		m_beliefs.clear();
		m_beliefs.resize(m_gmo.nvar(), factor(1.0));
		m_best_config.resize(m_gmo.nvar(), -1);

		// Output the bucket tree statistics
		std::cout << "Created bucket-tree with " << std::endl;
		std::cout << " - number of cliques:  " << C << std::endl;
		std::cout << " - number of edges:    " << elist.size() << std::endl;
		std::cout << " - max clique size:    " << max_clique_size << std::endl;
		std::cout << " - max separator size: " << max_sep_size << std::endl;
		std::cout << std::endl;

		if (m_debug) {
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
	/// \brief Compute the belief of a cluster.
	/// \param a 	The index of the cluster
	/// \return the factor representing the belief of the cluster.
	///
	factor calc_belief(findex a) {

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

	///
	/// \brief Compute the belief of a cluster excluding an incoming message.
	/// \param a 	The index of the cluster to compute the belief of
	/// \param i 	The index of the cluster sending the incoming message
	/// \return the factor representing the belief of cluster *a* excluding
	/// 	the incoming message from *i* to *a*.
	///
	factor incoming(findex a, size_t i) {

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

	///
	/// \brief Compute the belief of a cluster excluding backward messages.
	/// \param a 	The index of the cluster to compute the belief of
	/// \return the factor representing the belief of cluster *a* excluding
	/// 	the backward messages from clusters below *a*.
	///
	factor incoming(findex a) {

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

	///
	/// \brief Forward (top-down) message passing along the edge of the bucket tree.
	///
	void forward() {

		std::cout << "Begin forward (top-down) pass ..." << std::endl;

		m_log_z = 0; // reset the log parition function
		for (variable_order_t::const_iterator x = m_order.begin(); x != m_order.end(); ++x) {

			std::cout << " - Eliminating var " << *x
					<< (m_var_types[*x] ? " (MAP)\n" : " (SUM)\n");

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
		m_log_z += F.max();

		std::cout << "Finished forward pass with logZ: " << m_log_z << std::endl;
	}

	///
	/// \brief Backward (bottom-up) message passing along the edges of the bucket tree.
	///
	void backward() {

		// Calculate 'backward' messages
		std::cout << "Begin backward (bottom-up) pass ..." << std::endl;

		vector<std::pair<findex, findex> >::reverse_iterator ri = m_schedule.rbegin();
		for (; ri != m_schedule.rend(); ++ri ) {

			// Select backward message m(b->a)
			findex a = (*ri).first;
			findex b = (*ri).second;
			size_t i = m_edge_indeces[a][b]; // edge index

			// Select variables to eliminate
			variable_set VX = m_scopes[b] - m_separators[a][b];

			std::cout << " - Sending backward msg from " << a << " to " << b << std::endl;

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

		std::cout << "Finished backward pass." << std::endl;
	}

	///
	/// \brief Propagate the messages along the edges of the bucket tree.
	///
	void propagate() {

		std::cout << "Begin message passing over bucket tree ..." << std::endl;

		// Top-down 'forward' message passing
		forward();

		if (m_task == Task::PR) return; // stop here if we're only computing
										// the log partition function

		// Bottom-up 'backward' message passing
		backward();

		// Update beliefs/marginals
		update();
	}

	///
	/// \brief Update the beliefs (marginals or max-marginals)
	///
	void update() {

		// update beliefs (marginals) or compute the MAP/MMAP assignment
		if (m_task == Task::MAR) {
			for (vindex v = 0; v < m_gmo.nvar(); ++v) {
				findex c = m_clusters[v][0]; // get a cluster corresp. to current variable
				variable_set vars = m_scopes[c];
				variable VX = m_gmo.var(v);
				variable_set out = vars - VX;

				factor bel = calc_belief(c);
				m_beliefs[v] = marg(bel, VX);
				m_beliefs[v] /= std::exp(m_log_z); // normalize by logZ
			}
		} else if (m_task == Task::MAP) {
			for (variable_order_t::const_reverse_iterator x = m_order.rbegin();
					x != m_order.rend(); ++x) {

				variable VX = var(*x);
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
				variable VX = var(*x);
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


protected:
	// Members:

	graphical_model m_gmo; 				///< Original graphical model
	Task m_task;						///< Inference task
	OrderMethod m_order_method;			///< Variable ordering method
	size_t m_ibound;					///< Mini-bucket i-bound
	double m_log_z;						///< Log partition function value
	variable_order_t m_order;			///< Variable order
	std::vector<vindex> m_parents;		///< Pseudo tree
	vector<bool> m_var_types; 			///< Variable types (true if MAX, false if SUM)
	vector<factor> m_beliefs; 			///< Marginals
	std::vector<vindex> m_best_config;	///< MAP assignment
	std::vector<vindex> m_query; 		///< MAX variables for the MMAP task

private:
	// BT local structures:

	vector<bool> m_types;				///< The type of each cluster (SUM=false or MAX=true)
	vector<flist> m_clusters;			///< Clusters corresponding to each of the variables
	vector<flist> m_originals;			///< Original factors (index) for each cluster
	vector<variable_set> m_scopes;		///< The scope (vars) for each cluster
	vector<flist> m_in;					///< Incoming messages (source clusters) to each cluster
	vector<flist> m_out; 				///< Outgoing messages (target clusters) from each cluster
	flist m_roots;						///< Root cluster(s)
	vector<factor> m_forward; 			///< Forward messages (indexed by edge id in the bucket tree)
	vector<factor> m_backward; 			///< Backward messages (indexed by edge id in the bucket tree)

	vector<std::pair<findex, findex> > m_schedule;	///< Propagation schedule
	vector<vector<size_t> > m_edge_indeces;			///< Edge indeces
	vector<vector<variable_set> > m_separators; 	///< Separators between clusters

	std::map<size_t, size_t> m_cluster2var;			///< Maps cluster id to a variable id

	bool m_debug;						///< Internal flag for debugging only
};

} // namespace




#endif /* IBM_MERLIN_BE_H_ */
