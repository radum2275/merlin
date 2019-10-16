/*
 * wmb.h
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

/// \file wmb.h
/// \brief Weighted Mini-Buckets algorithm
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_MERLIN_WMB_H_
#define IBM_MERLIN_WMB_H_

#include "graphical_model.h"
#include "algorithm.h"

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

	~wmb() {};

	///
	/// \brief Clone the algorithm.
	/// \return the pointer to the new object containing the cloned algorithm.
	///
//	inline virtual wmb* clone() const {
//		wmb* gm = new wmb(*this);
//		return gm;
//	}

	// Can be an optimization algorithm or a summation algorithm....
	inline double ub() const {
		return m_logz;
	}
	inline double lb() const {
		throw std::runtime_error("Not implemented");
	}
	inline std::vector<size_t> best_config() const {
		return m_best_config;
	}

	inline double logZ() const {
		return m_logz;
	}
	inline double logZub() const {
		return m_logz;
	}
	inline double logZlb() const {
		return m_logz;
	}

	// No beliefs defined currently
	inline const factor& belief(size_t f) const {
		return m_beliefs[f];
	}
	inline const factor& belief(variable v) const {
		return m_beliefs[v];
	}
	inline const factor& belief(variable_set vs) const {
		throw std::runtime_error("Not implemented");
	}
	inline const vector<factor>& beliefs() const {
		return m_beliefs;
	}

	///
	/// \brief Get the original graphical model.
	///
	inline const graphical_model& get_gm_orig() const {
		return m_gmo;
	}

	///
	/// \brief Write the solution to the output stream.
	/// \param out		 		The output stream
	/// \param evidence 		The evidence variable value pairs
	/// \param old2new			The mapping between old and new variable indexing
	/// \param orig 			The graphical model prior to asserting evidence
	/// \param output_format	The output format (json or uai)
	///
	void write_solution(std::ostream& out, const std::map<size_t, size_t>& evidence,
			const std::map<size_t, size_t>& old2new, const graphical_model& orig,
			const std::set<size_t>& dummies, int output_format);

	///
	/// \brief Run the weighted mini-buckets algorithm.
	///
	void run();

	///
	/// \brief Inference tasks supported.
	///
	MER_ENUM( Task, PR,MAR,MAP,MMAP );

	///
	/// \brief Properties of the algorithm
	///
	MER_ENUM( Property , iBound,Order,Task,Iter,Debug,OrderIter );


	// Setting properties (directly or through property string):

	///
	/// \brief Set the mini-bucket i-bound parameter.
	///
	inline void set_ibound(size_t i) {
		m_ibound = i ? i : std::numeric_limits<size_t>::max();
	}

	///
	/// \brief Get the mini-bucket i-bound parameter.
	///
	inline size_t get_ibound() const {
		return m_ibound;
	}

	///
	/// \brief Set the variable types.
	///
	inline void set_var_types(const vector<bool>& var_types) {
		m_var_types = var_types;
	}

	///
	/// \brief Get the variable types.
	///
	inline const vector<bool>& get_var_types() const {
		return m_var_types;
	}

	/// 
	/// \brief Set the variable order.
	///
	inline void set_order(const variable_order_t& ord) {
		m_order = ord;
	}

	///
	/// \brief Set the variable order method.
	///
	inline void set_order_method(graphical_model::OrderMethod method) {
		m_order.clear();
		m_order_method = method;
	}

	///
	/// \brief Get the variable order.
	///
	inline const variable_order_t& get_order() {
		return m_order;
	}

	///
	/// \brief Get the pseudo tree.
	///
	inline const std::vector<vindex>& get_pseudo_tree() {
		return m_parents;
	}

	///
	/// \brief Set the pseudo tree.
	///
	inline void set_pseudo_tree(const vector<vindex>& p) {
		m_parents = p;
	}

	///
	/// \brief Set the query (MAP) variables.
	///
	inline void set_query(const std::vector<vindex>& q) {
		m_query = q;
	}

	///
	/// \brief Get the query (MAP) variables.
	///
	inline const std::vector<vindex>& get_query() {
		return m_query;
	}

	///
	/// \brief Set the graphical model.
	///
	inline void set_graphical_model(const graphical_model& gm) {
		m_gmo = gm;
	}

	///
	/// \brief Set the graphical model from a list of factors.
	///
	inline void set_graphical_model(const vector<factor>& fs) {
		m_gmo = graphical_model(fs);
	}

	///
	/// \brief Set the properties of the algorithm.
	/// \param opt 	The string containing comma separated property value pairs
	///	
	virtual void set_properties(std::string opt = std::string()) {
		if (opt.length() == 0) {
			set_properties("iBound=4,Order=MinFill,Iter=10,Task=MMAP,Debug=0,OrderIter=1");
			return;
		}
		m_debug = false;
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
			case Property::OrderIter:
				m_order_iter = atol(asgn[1].c_str());
				break;
			case Property::Debug:
				if (atol(asgn[1].c_str()) == 0) m_debug = false;
				else m_debug = true;
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
	inline factor elim(const factor& F, const variable_set& vs, const double w) {
		return F.sum_power(vs, w);
	}

	///
	/// \brief Compute the weighted marginal over a set of variables.
	/// \param F 	The reference of the factor to marginalize over
	/// \param vs 	The set of variables representing the scope of the marginal
	/// \param w 	The weight of the weighted elimination operator
	/// \return the factor representing the weighted marginal over the set of variables.
	///
	inline factor marg(const factor& F, const variable_set& vs, const double w) {
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
	inline double score(const vector<variable_set>& fin, const variable& VX, size_t i, size_t j) {
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
	void init();

	///
	/// \brief Compute the belief of a cluster.
	/// \param a 	The index of the cluster
	/// \return the factor representing the belief of the cluster.
	///
	factor calc_belief(findex a);

	///
	/// \brief Compute the belief of a cluster excluding an incoming message.
	/// \param a 	The index of the cluster to compute the belief of
	/// \param i 	The index of the cluster sending the incoming message
	/// \return the factor representing the belief of cluster *a* excluding
	/// 	the incoming message from *i* to *a*.
	///
	factor incoming(findex a, size_t i);

	///
	/// \brief Compute the belief of a cluster excluding backward messages.
	/// \param a 	The index of the cluster to compute the belief of
	/// \return the factor representing the belief of cluster *a* excluding
	/// 	the backward messages from clusters below *a*.
	///
	factor incoming(findex a);

	///
	/// \brief Forward (top-down) message passing with moment matching between the clusters of a bucket.
	///	
	void forward(double step);

	///
	/// \brief Backward (bottom-up) message passing.
	///
	void backward(size_t iter);

	///
	/// \brief Perform moment-matching between the clusters of a variable.
	///
	void match_clusters(size_t x, double step);

	///
	/// \brief Iterative tightening of the upper bound.
	///
	void tighten(size_t nIter, double stopTime = -1, double stopObj = -1);

	///
	/// \brief Update the beliefs (marginals or max-marginals)
	///
	void update();

protected:
	// Members:

	graphical_model m_gmo; 				///< Original graphical model
	Task m_task;						///< Inference task
	OrderMethod m_order_method;			///< Variable ordering method
	size_t m_order_iter;				///< Number of iterations for ordering heuristic
	size_t m_ibound;					///< Mini-bucket i-bound
	double m_logz;						///< Log partition function value
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
	std::map<size_t, size_t> m_cluster2var;			///< Maps cluster id to a variable id

	bool m_debug;						///< Internal debugging flag

};

} // namespace


#endif /* IBM_MERLIN_WMB_H_ */
