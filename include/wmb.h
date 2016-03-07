/*
 * wmb.h
 *
 *  Created on: 24 Mar 2015
 *      Author: radu
 */

/// \file wmb.h
/// \brief Weighted Mini-Buckets algorithm
/// \author Radu Marinescu

#ifndef IBM_MERLIN_WMB_H_
#define IBM_MERLIN_WMB_H_

#include "graphical_model.h"

namespace merlin {

/**
 * Weighted Mini-Buckets (WMB)
 *
 * Tasks supported: PR, MAR, MAP, MMAP
 *
 * WMB generalizes the classical MBE by replacing the sum operator with a
 * weighted sum operator, and thus it uses Holder's inequality to derive an
 * upper bound on the log partition function, MAP or Marginal MAP value. WMB
 * also uses an iterative cost-shifting scheme that matches marginals (weighted
 * or max) in order to tighten the upper-bound. Tightening is not guaranteed in
 * general, but it typically happens in practice.
 *
 */
class wmb: public graphical_model, public algorithm {
public:
	typedef graphical_model::findex findex;        ///< Factor index
	typedef graphical_model::vindex vindex;        ///< Variable index
	typedef graphical_model::flist flist;          ///< Collection of factor indices

public:

	///
	/// \brief Default constructor.
	///
	wmb() : graphical_model() {
		set_properties();
	}

	///
	/// \brief Constructor with a graphical model.
	///
	wmb(const graphical_model& gm) : graphical_model(gm), m_gmo(gm) {
		clear_factors();
		set_properties();
	}

	///
	/// \brief Clone the algorithm.
	/// \return the pointer to the new object containing the cloned algorithm.
	///
	virtual wmb* clone() const {
		wmb* gm = new wmb(*this);
		return gm;
	}

	// Can be an optimization algorithm or a summation algorithm....
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

	// No beliefs defined currently
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
	void write(const char* filename, const std::map<size_t, size_t>& evidence,
			const std::map<size_t, size_t>& old2new, const graphical_model& orig ) {

		// Open the output file
		std::ofstream out(filename);
		if (out.fail()) {
			throw std::runtime_error("Error while opening the output file.");
		}

		switch (m_task) {
		case Task::PR:
		case Task::MAR:
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
		case Task::MAP:
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
		case Task::MMAP:
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

	///
	/// \brief Run the weighted mini-buckets algorithm.
	///
	virtual void run() {
		init();
		tighten(m_num_iter);

		// Output solution (UAI output format)
		std::cout << "Converged after " << m_num_iter << " iterations in "
				<< (timeSystem() - m_start_time) << " seconds" << std::endl;

		switch (m_task) {
		case Task::PR:
		case Task::MAR:
			std::cout << "PR" << std::endl;
			std::cout << std::fixed << std::setprecision(6)
				<<m_log_z << " (" << std::scientific << std::setprecision(6)
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
		case Task::MAP:
			for (vindex v = 0; v < m_gmo.nvar(); ++v) {
				m_best_config[v] = m_beliefs[v].argmax();
			}
			m_lb = m_gmo.logP(m_best_config);
			std::cout << "Final Upper Bound is " << std::fixed << std::setw(12) << std::setprecision(6)
				<< m_log_z << " (" << std::scientific << std::setprecision(6)
				<< std::exp(m_log_z) << ")" << std::endl;
			std::cout << "Final Lower Bound is " << std::fixed << std::setw(12) << std::setprecision(6)
				<< m_lb << " (" << std::scientific << std::setprecision(6)
				<< std::exp(m_lb) << ")" << std::endl;
			std::cout << "MAP" << std::endl;
			std::cout << m_gmo.nvar();
			for (vindex v = 0; v < m_gmo.nvar(); ++v) {
				//std::cout << " " << m_beliefs[v].argmax();
				std::cout << " " << m_best_config[v];
			}
			std::cout << std::endl;

			break;
		case Task::MMAP:
			for (size_t i = 0; i < m_query.size(); ++i) {
				vindex v = m_query[i];
				m_best_config[v] = m_beliefs[v].argmax();
			}
			std::cout << "Final Upper Bound is " << std::fixed << std::setw(12) << std::setprecision(6)
				<< m_log_z << " (" << std::scientific << std::setprecision(6)
				<< std::exp(m_log_z) << ")" << std::endl;
			std::cout << "MMAP" << std::endl;
			std::cout << m_query.size();
			for (vindex v = 0; v < m_gmo.nvar(); ++v) {
				if (m_var_types[v] == true)
					std::cout << " " << m_best_config[v];
			}
			std::cout << std::endl;
		}

	}

	///
	/// \brief Inference tasks supported.
	///
	MER_ENUM( Task, PR,MAR,MAP,MMAP );

	///
	/// \brief Properties of the algorithm
	///
	MER_ENUM( Property , iBound,Order,Task,Iter );


	// Setting properties (directly or through property string):

	///
	/// \brief Set the mini-bucket i-bound parameter.
	///
	void set_ibound(size_t i) {
		m_ibound = i ? i : std::numeric_limits<size_t>::max();
	}

	///
	/// \brief Get the mini-bucket i-bound parameter.
	///
	size_t get_ibound() const {
		return m_ibound;
	}

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
			set_properties("iBound=4,Order=MinFill,Iter=100,Task=MMAP");
			return;
		}
		std::vector<std::string> strs = merlin::split(opt, ',');
		for (size_t i = 0; i < strs.size(); ++i) {
			std::vector<std::string> asgn = merlin::split(strs[i], '=');
			switch (Property(asgn[0].c_str())) {
			case Property::iBound:
				set_ibound(atol(asgn[1].c_str()));
				break;
			case Property::Order:
				m_order.clear();
				m_parents.clear();
				m_order_method = graphical_model::OrderMethod(asgn[1].c_str());
				break;
			case Property::Task:
				m_task = Task(asgn[1].c_str());
				break;
			case Property::Iter:
				m_num_iter = atol(asgn[1].c_str());
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
	factor marg(const factor& F, const variable_set& vs, const double w) {
		return F.marginal(vs, w);
	}

	///
	/// \brief Scoring function for bucket aggregation.
	/// \param fin 		The set of factor scopes containing 
	///						the pair (i,j) to be aggregated
	/// \param VX 		The bucket variable
	/// \param i 		The index of first scope
	/// \param j 		The index of the second pair
	/// \return the score that corresponds to aggregating the two scopes.
	///		It returns -3 if unable to combine, -1 for scope only aggregation,
	///		and otherwise a positive double score.
	///
	double score(const vector<variable_set>& fin, const variable& VX, size_t i, size_t j) {
		double err;
		const variable_set& F1 = fin[i], &F2 = fin[j];           // (useful shorthand)
		size_t iBound = std::max(std::max(m_ibound, F1.nvar() - 1),
				F2.nvar() - 1);      // always OK to keep same size
		variable_set both = F1 + F2;
		if (both.nvar() > iBound+1)
			err = -3;  // too large => -3
		else
			err = 1.0 / (F1.nvar() + F2.nvar()); // greedy scope-based 2 (check if useful???)
		//else if (_byScope) err = 1;            // scope-based => constant score
		return err;
	}

	///
	/// \brief Helper class for pairs of sorted indices.
	///
	struct Pair: public std::pair<size_t, size_t> {
		Pair(size_t ii, size_t jj) {
			if (ii < jj) {
				first = jj;
				second = ii;
			} else {
				first = ii;
				second = jj;
			}
		}
	};

	///
	/// \brief Initialize the weighted mini-buckets algorithm.
	///
	void init() {

		// Start the timer and store it
		m_start_time = timeSystem();

		// Initialize variable types
		m_var_types.resize(m_gmo.nvar(), false); // all SUM
		for (size_t i = 0; i < m_query.size(); ++i) {
			m_var_types[m_query[i]] = true; // set MAP variable
		}

		// Prologue
		std::cout << VERSIONINFO << std::endl << COPYRIGHT << std::endl;
		std::cout << "Initialize inference engine ..." << std::endl;
		std::cout << "+ tasks supported  : PR, MAR, MAP, MMAP" << std::endl;
		std::cout << "+ algorithm        : " << "WMB" << std::endl;
		std::cout << "+ i-bound          : " << m_ibound << std::endl;
		std::cout << "+ iterations       : " << m_num_iter << std::endl;
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
			m_parents.clear(); // (new elim order => need new pseudotree) !!! should do together
			std::copy(m_order.begin(), m_order.end(),	std::ostream_iterator<size_t>(std::cout, " "));
		}
		if (m_parents.size() == 0) {     // if we need to construct a pseudo-tree
			m_parents = m_gmo.pseudo_tree(m_order);
		}

		std::cout << std::endl;
		size_t wstar = m_gmo.induced_width(m_order);
		std::cout << "+ induced width    : " << wstar << std::endl;
		std::cout << "+ exact inference  : " << (m_ibound >= wstar ? "Yes" : "No") << std::endl;
		if (m_ibound >= wstar) m_num_iter = 1; // exact inference requires 1 iteration over the join-tree

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
			double weight = 1.0/((double)ids.size()); // uniform weights
			// check if bucket variable is a MAP variable
			if (m_var_types[*x] == true) weight = infty();

			// Eliminate individually each mini-bucket
			vector<findex> alphas;
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

				// now incoming nodes to *i is just alpha
				Orig[*i].clear();
				New[*i].clear();
				New[*i] |= alpha;

				// recompute and update adjacency
				insert(vin, *i, fin[*i]);
			}

			//std::cout<<"\n";

		}
		// end for: variable elim order

		// separators and cluster scopes
		size_t C = m_factors.size();
		m_separators.resize(C);
		for (size_t i = 0; i < C; ++i) m_separators[i].resize(C);
		m_scopes.resize(C);
		for (size_t i = 0; i < C; ++i) m_scopes[i] = m_factors[i].vars();
		const std::vector<edge_id>& elist = edges();
		for (size_t i = 0; i < elist.size(); ++i) {
			findex a,b;
			a = elist[i].first;
			b = elist[i].second;
			if (a > b) continue;
			variable_set sep = m_factors[a].vars() & m_factors[b].vars();
			m_separators[a][b] = sep;
			m_separators[b][a] = sep;
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
			m_factors[i] = get_factor(1.0); // init

			// clique potential
			for (flist::const_iterator j = m_originals[i].begin();
					j != m_originals[i].end(); ++j) {
				m_factors[i] *= m_gmo.get_factor(*j);
			}

		}

		// initialize beliefs (marginals)
		m_log_z = 0;
		m_beliefs.clear();
		m_beliefs.resize(m_gmo.nvar(), get_factor(1.0));
		m_reparam.resize( m_factors.size(), get_factor(1.0) );
		m_best_config.resize(m_gmo.nvar(), -1);
	}

	///
	/// \brief Compute the belief of a cluster.
	/// \param a 	The index of the cluster
	/// \return the factor representing the belief of the cluster.
	///
	factor calc_belief(findex a) {

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

	///
	/// \brief Compute the belief of a cluster excluding an incoming message.
	/// \param a 	The index of the cluster to compute the belief of
	/// \param i 	The index of the cluster sending the incoming message
	/// \return the factor representing the belief of cluster *a* excluding
	/// 	the incoming message from *i* to *a*.
	///
	factor incoming(findex a, size_t i) {

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

	///
	/// \brief Compute the belief of a cluster excluding backward messages.
	/// \param a 	The index of the cluster to compute the belief of
	/// \return the factor representing the belief of cluster *a* excluding
	/// 	the backward messages from clusters below *a*.
	///
	factor incoming(findex a) {

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

	// forward pass with moment matching (between the clusters of a bucket)
	///
	/// \brief Forward (top-down) message passing with moment matching between the clusters of a bucket.
	///	
	void forward(double step) {

		m_log_z = 0;
		for (variable_order_t::const_iterator x = m_order.begin(); x != m_order.end(); ++x) {

//			std::cout << "Eliminating "<<*x << (_varTypes[*x] ? "(MAP)\n" : "(SUM)\n");

			// moment-match the clusters of this bucket
			match_clusters(*x, step);

			// generate forward messages
			variable VX = var(*x);
			for (flist::const_iterator it = m_clusters[*x].begin();
					it != m_clusters[*x].end(); ++it) {
				findex a = (*it);
				if ( m_out[a].size() > 0 ) {
					findex b = *(m_out[a].begin());
					size_t i = m_edge_indeces[a][b];

//					std::cout << "forward msg (" << a << "," << b << "): elim = " << VX << " -> ";

					factor tmp = incoming(a, i);
					if (m_var_types[*x] == false) { // SUM
						m_forward[i] = tmp.sum_power(VX, 1.0/m_weights[a]);
						m_forward[i] /= m_forward[i].sum();
					} else { // MAX
						m_forward[i] = tmp.max(VX);
						double mx = m_forward[i].max();
						m_forward[i] /= mx;
						m_log_z += std::log(mx);
					}

//					std::cout << _forward[i] << std::endl;
				}
			}
		}

		// compute the upper bound
		factor F(0.0);
		for (flist::const_iterator ci = m_roots.begin();
				ci != m_roots.end(); ++ci) {

			factor bel = calc_belief(*ci);
			if (m_task == Task::MAP)
				F += log( bel.max());
			else
				F += log( bel.sum() );
		}
		m_log_z += F.max();
	}

	///
	/// \brief Backward (bottom-up) message passing.
	///
	void backward(size_t iter) {

		// update backward messages
		vector<std::pair<findex, findex> >::reverse_iterator ri = m_schedule.rbegin();
		for (; ri != m_schedule.rend(); ++ri ) {

			// compute backward message m(b->a)
			findex a = (*ri).first;
			findex b = (*ri).second;
			size_t i = m_edge_indeces[a][b]; // edge index

			variable_set VX = m_scopes[b] - m_separators[a][b];

//			std::cout << "backward msg (" << b << "," << a << "): elim = " << VX << " -> ";

			// compute the belief at b
			factor bel = calc_belief(b);

			if (m_types[b] == false && m_types[a] == false) { // SUM-SUM

				bel ^= 1.0/m_weights[b];
				bel /= (m_forward[i]^(1.0/m_weights[a])); // divide out m(a->b)

				//_backward[i] = elim(bel, VX, 1);
				m_backward[i] = bel.sum(VX);
				m_backward[i] ^= (m_weights[b]);
				m_backward[i] /= m_backward[i].sum(); // normalize

			} else if (m_types[b] == true && m_types[a] == true) { // MAX-MAX

				bel /= m_forward[i]; // divide out m(a->b)
				m_backward[i] = bel.max(VX);
				m_backward[i] /= m_backward[i].max(); // normalize

			} else if (m_types[b] == true && m_types[a] == false) { // MAX-SUM

				bel = bel.sigma(iter); // the sigma operator that focuses on max
				bel /= (m_forward[i]^(1.0/m_weights[a])); // divide out m(a->b)

				m_backward[i] = bel.sum(VX);
				m_backward[i] ^= (m_weights[a]);
				m_backward[i] /= m_backward[i].max(); // normalize

			} else {
				assert(false); // cannot reach this case!!
			}

//			std::cout << _backward[i] << std::endl;

			// normalize
//			double mx = m_backward[i].max();
//			_backwardMsg[i] /= mx;
		}

	}

	///
	/// \brief Perform moment-matching between the clusters of a variable.
	///
	void match_clusters(size_t x, double step) {

		if (m_clusters[x].size() <= 1)
			return; // no matching

		variable VX = var(x);
		if (m_var_types[x] == true) { // max marginals matching

//			std::cout << "matching max marginals" << std::endl;

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

//				std::cout << " reparam      : " << m_reparam[a] << std::endl;
			}

		} else { // weighted marginals matching

//			std::cout << "matching weighted marginals on cliques: ";
//			std::copy(_match[x].begin(), _match[x].end(),
//					std::ostream_iterator<int>(std::cout, " "));
//			std::cout << std::endl;

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

//				std::cout << " clique belief: " << bel << std::endl;
//				std::cout << " marginal     : " << ftmp[i] << std::endl;
			}

//			std::cout << " geom mean    : " << fmatch << std::endl;
			i = 0;
			for (flist::const_iterator it = m_clusters[x].begin();
					it != m_clusters[x].end(); ++it, ++i) {
				findex a = (*it);
				m_reparam[a] *= ((fmatch/ftmp[i])^(step*m_weights[a]));

//				std::cout << " reparam      : " << m_reparam[a] << std::endl;
			}
		}
	}

	///
	/// \brief Iterative tightering of the upper bound.
	///
	void tighten(size_t nIter, double stopTime = -1, double stopObj = -1) {
		std::cout << "Begin message passing over join graph ..." << std::endl;
		std::cout << "- stopObj  : " << stopObj << std::endl;
		std::cout << "- stopTime : " << stopTime << std::endl;
		std::cout << "- stopIter : " << nIter << std::endl;

		double minZ = infty();
		vector<factor> minB;
		for (size_t iter = 1; iter <= nIter; ++iter) {
			double step = 1.0/(double)iter;
			double prevZ = m_log_z;

			forward(step);
			backward(iter);
			update();

			// keep track of tightest upper bound
			if (m_log_z < minZ) {
				minZ = m_log_z;
				minB = m_beliefs;
			}

			double dObj = fabs(m_log_z - prevZ);
			std::cout << "  WMB: " << std::fixed << std::setw(12) << std::setprecision(6)
				<< m_log_z << " (" << std::scientific << std::setprecision(6)
				<< std::exp(m_log_z) << ") ";
			std::cout << "\td=" << dObj << "\t time="  << std::fixed
				<< std::setprecision(6) << (timeSystem() - m_start_time)
				<< "\ti=" << iter << std::endl;

			if (dObj < stopObj) break;

			// do at least one iterations
			if (stopTime > 0 && stopTime <= (timeSystem() - m_start_time))
				break;
		} // end for

		m_log_z = minZ; // keep tightest upper bound
		m_beliefs = minB; // and corresponding beliefs (marginals)
	}

	///
	/// \brief Update the beliefs (marginals or max-marginals)
	///
	void update() {

		// update beliefs (marginals)
		for (vindex v = 0; v < m_gmo.nvar(); ++v) {
			findex c = m_clusters[v][0]; // get a cluster corresp. to current variable
			double w = m_weights[c];
			variable_set vars = m_scopes[c];
			variable VX = m_gmo.var(v);
			variable_set out = vars - VX;

			factor bel = calc_belief(c);
			m_beliefs[v] = marg(bel, VX, w);
			m_beliefs[v] /= m_beliefs[v].sum(); // normalize
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
	size_t m_num_iter; 					///< Number of iterations to be executed
	double m_lb;						///< Lower bound (ie, value of MAP assignment)

private:
	// JG local structures:

	vector<bool> m_types;				///< The type of each cluster (SUM=false or MAX=true)
	vector<double> m_weights;			///< The weight of each cluster
	vector<flist> m_clusters;			///< Clusters to be matched for each variable
	vector<flist> m_originals;			///< Original factors (index) for each cluster
	vector<variable_set> m_scopes;		///< The scope (vars) for each cluster
	vector<flist> m_in;					///< Incoming to each cluster
	vector<flist> m_out; 				///< Outgoing from each cluster
	flist m_roots;						///< Root cluster(s)
	vector<factor> m_forward; 			///< Forward messages (by edge)
	vector<factor> m_backward; 			///< Backward messages (by edge)
	vector<factor> m_reparam; 			///< Reparameterization function (by cluster)

	vector<std::pair<findex, findex> > m_schedule;	///< Propagation schedule
	vector<vector<size_t> > m_edge_indeces;			///< Edge indeces
	vector<vector<variable_set> > m_separators; 	///< Separators between clusters

};

} // namespace


#endif /* LIB_IJGP_H_ */
