/*
 * ijgp.h
 *
 *  Created on: 24 Mar 2015
 *      Author: radu
 */

/// \file ijgp.h
/// \brief Iterative Join Graph Propagation (IJGP) algorithm
/// \author Radu Marinescu

#ifndef IBM_MERLIN_IJGP_H_
#define IBM_MERLIN_IJGP_H_

#include "graphical_model.h"

namespace merlin {

/**
 * Iterative Join-Graph Propagation (IJGP)
 *
 * Based on [Dechter and Mateescu, 2002] and [Marinescu, Kask and Dechter, 2003]
 *
 * Tasks supported: MAR and MAP
 *
 * IJGP is parameterized by an i-bound which limits the size of each cluster in
 * the join-graph to at most i distict variables. Clearly IJGP(1) is
 * equivalent with Loopy Belief Propagation, while IJGP(w*) is equivalent with
 * the Join-Tree algorithm, hence exact.
 *
 * The join-graph used by IJGP is obtained by running the mini-bucket algorithm
 * schematically (ie, without computing the actual messages, only their scopes)
 * and then connecting the mini-buckets residing in the same bucket. Messages
 * are then propagated along the join-graph edges, following a top-down or
 * bottom-up schedule.
 *
 * Note that IJGP is only used for MAR (sum-prod) and MAP (max-prod) tasks. It
 * doesn't compute an upper-bound (on the partition function, or the MAP value)
 * because of overcounting. Therefore, logZ reported during the execution of
 * the algorithm shouldn't be used as a valid measure for bounding (ignore). For
 * valid bounding, use the WMB algorithm implemented in this library.
 */

class ijgp: public graphical_model, public algorithm {
public:
	typedef graphical_model::findex findex;      ///< Factor index
	typedef graphical_model::vindex vindex;      ///< Variable index
	typedef graphical_model::flist flist;		///< Collection of factor indices

public:

	///
	/// \brief Default constructor.
	///
	ijgp() : graphical_model() {
		set_properties();
	}
	
	///
	/// \brief Constructor.
	///
	ijgp(const graphical_model& gm) : graphical_model(gm), m_gmo(gm) {
		clear_factors();
		set_properties();
	}

	///
	/// \brief Clone the model.
	/// \return the pointer to the object containing the cloned model.
	///
	virtual ijgp* clone() const {
		ijgp* gm = new ijgp(*this);
		return gm;
	}

	// Can be an optimization algorithm or a summation algorithm

	double ub() const {
		throw std::runtime_error("IJGP does not compute an upper bound due to overcounting.");
	}
	double lb() const {
		throw std::runtime_error("IJGP does not compute a lower bound due to overcounting.");
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
	const std::vector<factor>& beliefs() const {
		return m_beliefs;
	}

	///
	/// \brief Access the original graphical model.
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
		case Task::MAR:
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
		}
	}

	virtual void run() {

		init();
		propagate(num_iter);

		// Output solution (UAI output format)
		std::cout << "Converged after " << num_iter << " iterations in "
				<< (timeSystem() - m_start_time) << " seconds" << std::endl;

		switch (m_task) {
		case Task::MAR:
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
			m_best_config.resize(m_gmo.nvar());
			for (vindex v = 0; v < m_gmo.nvar(); ++v) {
				m_best_config[v] = m_beliefs[v].argmax();
			}
			m_lb = m_gmo.logP(m_best_config);
			std::cout << "Final Lower Bound is " << std::fixed << std::setw(12) << std::setprecision(6)
				<< m_lb << " (" << std::scientific << std::setprecision(6)
				<< std::exp(m_lb) << ")" << std::endl;
			std::cout << "MAP" << std::endl;
			std::cout << m_gmo.nvar();
			for (vindex v = 0; v < m_gmo.nvar(); ++v) {
				std::cout << " " << m_beliefs[v].argmax();
			}
			std::cout << std::endl;

			break;
		}
	}

	///
	/// \brief Inference tasks supported.
	///
	MER_ENUM( Task, MAR,MAP );

	///
	/// \brief Properties of the algorithm.
	///
	MER_ENUM( Property , iBound,Order,Iter,Task );

	///
	/// \brief Elimination operators (sum, max).
	///
	MER_ENUM( ElimOp , Max,Sum );


protected:
	// Members:

	graphical_model m_gmo;					///< Original graphical model.
	size_t num_iter;						///< Number of iterations
	Task m_task;							///< Inference task
	ElimOp m_elim_op;						///< Elimination operator
	size_t m_ibound;						///< i-bound parameter
	double m_log_z;							///< Log partition function value
	variable_order_t m_order;				///< Variable elimination order
	OrderMethod m_order_method;				///< Ordering method
	std::vector<vindex> m_parents;			///< Pseudo tree
	std::vector<factor> m_beliefs; 			///< Marginals (or beliefs)
	std::vector<size_t> m_best_config;		///< MAP assignment
	double m_lb; 							///< Lower bound (ie, value of the MAP assignment)

private:

	// JG local structures:
	vector<flist> m_clusters;					///< Clusters for each variable
	vector<vector<variable_set> > m_separators; ///< Separators between clusters
	vector<flist> m_originals;					///< Original factors (index) for each cluster
	vector<variable_set> m_scopes;				///< The scope (vars) for each cluster
	vector<flist> m_in;							///< Incoming to each cluster
	vector<flist> m_out; 						///< Outgoing from each cluster
	flist m_roots;								///< Root cluster(s)
	vector<factor> m_forward;					///< Forward messages (by edge)
	vector<factor> m_backward; 					///< Backward messages (by edge)
	vector<std::pair<findex, findex> > m_schedule; ///< Propagation schedule
	vector<vector<size_t> > m_edge_indeces;		///< Edge indexes

public:
	// Setting properties (directly or through property string):

	///
	/// \brief Set the i-bound parameter.
	/// \param i 	The value of the i-bound (> 1)
	///
	void set_ibound(size_t i) {
		m_ibound = i ? i : std::numeric_limits<size_t>::max();
	}

	///
	/// \brief Return the i-bound parameter.
	///
	size_t get_ibound() const {
		return m_ibound;
	}

	///
	/// \brief Set the variable elimination order.
	/// \param ord 	The variable order
	///
	void set_order(const variable_order_t& ord) {
		m_order = ord;
	}

	///
	/// \brief Set the variable elimination order method.
	/// \param method 	The elimination order method
	///
	void set_order_method(OrderMethod method) {
		m_order.clear();
		m_order_method = method;
	}

	///
	/// \brief Return the variable elimination order.
	///
	const variable_order_t& get_order() {
		return m_order;
	}

	///
	/// \brief Return the pseudo tree.
	///
	const std::vector<vindex>& get_pseudo_tree() {
		return m_parents;
	}

	///
	/// \brief Set the pseudo tree.
	/// \param p 	The vector representing the pseudo tree
	///
	void set_pseudo_tree(const vector<vindex>& p) {
		m_parents = p;
	}

	///
	/// \brief Set the graphical model content.
	/// \param gm 	The reference to the graphical model object
	///
	void set_graphical_model(const graphical_model& gm) {
		m_gmo = gm;
	}

	///
	/// \brief Set the graphical model content from a list of factors.
	/// \param fs 	The list of factors
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
			set_properties("iBound=4,Order=MinFill,Iter=10,Task=MAR");
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
			case Property::Iter:
				num_iter = atol(asgn[1].c_str());
				break;
			case Property::Task:
				m_task = Task(asgn[1].c_str());
				if (m_task == Task::MAR) m_elim_op = ElimOp::Sum;
				else m_elim_op = ElimOp::Max;
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
	/// \return the factor resulted from eliminating the set of variables.
	///
	factor elim(const factor& F, const variable_set& vs) {
		switch ( m_elim_op ) {
		case ElimOp::Sum:
			return F.sum(vs);
			break;
		case ElimOp::Max:
			return F.max(vs);
			break;
		}
		throw std::runtime_error("Unknown elim op");
	}

	///
	/// \brief Compute the marginal over a set of variables.
	/// \param F 	The reference of the factor to marginalize over
	/// \param vs 	The set of variables representing the scope of the marginal
	/// \return the factor representing the marginal over the set of variables.
	///
	factor marg(const factor& F, const variable_set& vs) {
		switch ( m_elim_op ) {
		case ElimOp::Sum:
			return F.marginal(vs);
			break;
		case ElimOp::Max:
			return F.maxmarginal(vs);
			break;
		}
		throw std::runtime_error("Unknown elim op");
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
	struct sPair: public std::pair<size_t, size_t> {
		sPair(size_t ii, size_t jj) {
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
	/// \brief Create the mini-bucket based join-graph (symbols only).
	///
	void init() {

		// Start the timer and store it
		m_start_time = timeSystem();

		// Prologue
		std::cout << VERSIONINFO << std::endl << COPYRIGHT << std::endl;
		std::cout << "Initialize inference engine ..." << std::endl;
		std::cout << "+ tasks supported  : MAR, MAP" << std::endl;
		std::cout << "+ algorithm        : " << "IJGP" << std::endl;
		std::cout << "+ i-bound          : " << m_ibound << std::endl;
		std::cout << "+ iterations       : " << num_iter << std::endl;
		std::cout << "+ inference task   : " << m_task << std::endl;
		std::cout << "+ ordering heur.   : " << m_order_method << std::endl;
		std::cout << "+ elimination      : ";

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
		std::cout << "+ induced width    : " << wstar << std::endl;
		std::cout << "+ exact inference  : " << (m_ibound >= wstar ? "Yes" : "No") << std::endl;
		if (m_ibound >= wstar) num_iter = 1; // exact inference requires 1 iteration over the join-tree

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

			//std::cout << "Eliminating " << *x << std::endl;

			variable VX = var(*x);
			if (*x >= vin.size() || vin[*x].size() == 0)
				continue;  // check that we have some factors over this variable

			flist ids = vin[*x];  // list of factor IDs contained in this bucket

			//// Select allocation into mini-buckets ///////////////////////////
			typedef flist::const_iterator flistIt;
			typedef std::pair<double, sPair> _INS;
			std::multimap<double, sPair> scores;
			std::map<sPair, std::multimap<double, sPair>::iterator> reverseScore;

			//// Populate list of pairwise scores for aggregation //////////////
			for (flistIt i = ids.begin(); i != ids.end(); ++i) {
				for (flistIt j = ids.begin(); j != i; ++j) {
					double err = score(fin, VX, *i, *j);
					sPair sp(*i, *j);
					reverseScore[sp] = scores.insert(_INS(err, sp)); // save score
				}
				reverseScore[sPair(*i, *i)] = scores.insert(
						_INS(-1, sPair(*i, *i)));       // mark self index at -1
			}

			//// Run through until no more pairs can be aggregated: ////////////
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

			//std::cout << " - mini-buckets: " << ids.size() << std::endl;

			//// Eliminate individually each mini-bucket ///////////////////////
			vector<findex> alphas;
			for (flistIt i = ids.begin(); i != ids.end(); ++i) {
				//
				// Create new cluster alpha over this set of variables; save function parameters also
				findex alpha = findex(-1);
				alpha = add_factor(factor(fin[*i]));
				alphas.push_back(alpha);
				m_clusters[*x] |= alpha;

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
			//std::cout<<"\n";

		}
		/// end for: variable elim order ///////////////////////////////////////

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
		m_beliefs.resize(m_gmo.nvar(), factor(1.0));

		// output summary of initialization
		std::cout << "Created join graph with " << m_factors.size() << " clusters and "
				<< elist.size() << " edges" << std::endl;

	}

	///
	/// \brief Compute the belief of a cluster.
	/// \param a 	The index of the cluster
	/// \return the factor representing the belief of the cluster.
	///
	factor _belief(findex a) {

		factor bel = m_factors[a];

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
	/// \param b 	The index of the cluster sending the incoming message
	/// \return the factor representing the belief of cluster *a* excluding
	/// 	the incoming message from *b* to *a*.
	///
	factor _belief(findex a, size_t b) {

		factor bel = m_factors[a];

		// forward messages to 'a'
		for (flist::const_iterator ci = m_in[a].begin();
				ci != m_in[a].end(); ++ci) {
			findex p = (*ci);
			if (p == b) continue;
			size_t j = m_edge_indeces[p][a];
			bel *= m_forward[j];
		}

		// backward message to 'a'
		for (flist::const_iterator ci = m_out[a].begin();
				ci != m_out[a].end(); ++ci) {
			findex p = (*ci);
			if (p == b) continue;
			size_t j = m_edge_indeces[a][p];
			bel *= m_backward[j];
		}

		return bel;
	}


	///
	/// \brief Forward (top-down) message passing.
	///
	void forward() {

		// compute forward messages
		//std::cout << "FORWARD MESSAGES" << std::endl;
		vector<std::pair<findex, findex> >::iterator fi = m_schedule.begin();
		for (; fi != m_schedule.end(); ++fi ) {

			// compute forward message m(a->b)
			findex a = (*fi).first;  // cluster index a
			findex b = (*fi).second; // cluster index b
			size_t ei = m_edge_indeces[a][b]; // edge index

			// variables to eliminate
			variable_set VX = m_scopes[a] - m_separators[a][b];

			//std::cout << "  forward msg (" << a << "," << b << "): elim = " << VX << " -> ";

			// compute the belief at a (excluding message b->a)
			factor bel = _belief(a, b);
			m_forward[ei] = elim(bel, VX); //bel.sum(VX);
			if (m_task == Task::MAP) {
				double mx = m_forward[ei].max();
				m_forward[ei] /= mx;
				m_log_z += std::log(mx);
			} else {
				m_forward[ei] /= m_forward[ei].sum();
			}

			//std::cout << _forwardMsg[ei] << std::endl;
		}
	}

	///
	/// \brief Backward (bottom-up) message passing.
	///
	void backward() {

		// update backward messages
		//std::cout << "BACKWARD MESSAGES:" << std::endl;
		vector<std::pair<findex, findex> >::reverse_iterator ri = m_schedule.rbegin();
		for (; ri != m_schedule.rend(); ++ri ) {

			// compute backward message m(b->a)
			findex a = (*ri).first;
			findex b = (*ri).second;
			size_t ei = m_edge_indeces[a][b]; // edge index

			variable_set VX = m_scopes[b] - m_separators[a][b];

			//std::cout << "  backward msg (" << b << "," << a << "): elim = " << VX << " -> ";

			// compute the belief at b (excluding message a->b)
			factor bel = _belief(b, a);
			m_backward[ei] = elim(bel, VX); //bel.sum(VX);
			if (m_task == Task::MAP) {
				double mx = m_backward[ei].max();
				m_backward[ei] /= mx;
				m_log_z += std::log(mx);
			} else {
				m_backward[ei] /= m_backward[ei].sum();
			}

			//std::cout << _backwardMsg[ei] << std::endl;
		}

	}

	///
	/// \brief Update the beliefs (marginals or max-marginals) for each variable.
	///
	void update() {

		// update beliefs
		for (vindex v = 0; v < m_gmo.nvar(); ++v) {
			findex c = m_clusters[v][0]; // get a cluster corresp. to current variable
			variable_set vars = m_scopes[c];
			variable VX = m_gmo.var(v);
			variable_set out = vars - VX;

			factor bel = _belief(c);
			m_beliefs[v] = marg(bel, VX);
			if (m_task == Task::MAP) {
				double mx = m_beliefs[v].max();
				m_beliefs[v] /= mx;
				m_log_z += std::log(mx);
			} else {
				m_beliefs[v] /= m_beliefs[v].sum(); // normalize
			}
		}

		// compute logZ
		factor F(0.0);
		for (flist::const_iterator ci = m_roots.begin();
				ci != m_roots.end(); ++ci) {

			factor bel = _belief(*ci);
			if (m_task == Task::MAP) {
				F += log( bel.max() );
			} else {
				F += log( bel.sum() );
			}
		}
		m_log_z = F.max();
	}

	///
	/// \brief Iterative message passing over the join graph.
	/// \param nIter 	The number of iterations
	/// \param stopTime	The time limit
	/// \param stopObj 	The error tolerance (ie, difference between objective values)
	///
	void propagate(size_t nIter, double stopTime = -1, double stopObj = -1) {
		std::cout << "Begin message passing over join graph ..." << std::endl;
		std::cout << "- stopObj  : " << stopObj << std::endl;
		std::cout << "- stopTime : " << stopTime << std::endl;
		std::cout << "- stopIter : " << nIter << std::endl;

		for (size_t iter = 1; iter <= nIter; ++iter) {

			double prevZ = m_log_z;
			m_log_z = 0;
			forward();
			backward();
			update(); // update beliefs

			double dObj = fabs(m_log_z - prevZ);
			std::cout << "  IJGP: " << std::fixed << std::setw(12) << std::setprecision(6)
				<< m_log_z << " (" << std::scientific << std::setprecision(6)
				<< std::exp(m_log_z) << ") ";
			std::cout << "\td=" << dObj << "\t time="  << std::fixed
				<< std::setprecision(6) << (timeSystem() - m_start_time)
				<< "\ti=" << iter << std::endl;

			if (dObj < stopObj) break;

			// do at least one iterations
			if (stopTime > 0 && stopTime <= (timeSystem() - m_start_time))
				break;
		}
	}


};

} // namespace


#endif /* LIB_IJGP_H_ */
