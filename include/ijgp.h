/*
 * ijgp.h
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

/// \file ijgp.h
/// \brief Iterative Join Graph Propagation (IJGP) algorithm
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_MERLIN_IJGP_H_
#define IBM_MERLIN_IJGP_H_

#include "graphical_model.h"
#include "algorithm.h"

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
//	virtual ijgp* clone() const {
//		ijgp* gm = new ijgp(*this);
//		return gm;
//	}

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
		return m_logz;
	}
	double logZub() const {
		return m_logz;
	}
	double logZlb() const {
		return m_logz;
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
	/// \brief Run the algorithm
	///
	void run();

	///
	/// \brief Inference tasks supported.
	///
	MER_ENUM( Task, PR,MAR,MAP );

	///
	/// \brief Properties of the algorithm.
	///
	MER_ENUM( Property , iBound,Order,Iter,Task,Debug );

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
	double m_logz;							///< Log partition function value
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
	std::map<size_t, size_t> m_cluster2var;		///< Maps cluster id to a variable id

	bool m_debug;								///< Internal debugging flag

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
			set_properties("iBound=4,Order=MinFill,Iter=10,Task=MAR,Debug=0");
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
			case Property::Iter:
				num_iter = atol(asgn[1].c_str());
				break;
			case Property::Task:
				m_task = Task(asgn[1].c_str());
				if (m_task == Task::MAR) m_elim_op = ElimOp::Sum;
				else m_elim_op = ElimOp::Max;
				break;
			case Property::Debug:
				m_debug = (atol(asgn[1].c_str()) == 0) ? false : true;
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
	/// \param b 	The index of the cluster sending the incoming message
	/// \return the factor representing the belief of cluster *a* excluding
	/// 	the incoming message from *b* to *a*.
	///
	factor calc_belief(findex a, size_t b);

	///
	/// \brief Compute the belief of a cluster excluding backward messages.
	/// \param a 	The index of the cluster to compute the belief of
	/// \return the factor representing the belief of cluster *a* excluding
	/// 	the backward messages from clusters below *a*.
	///
	factor incoming(findex a);

	///
	/// \brief Forward (top-down) message passing.
	///
	void forward();

	///
	/// \brief Backward (bottom-up) message passing.
	///
	void backward();

	///
	/// \brief Update the beliefs (marginals or max-marginals) for each variable.
	///
	void update();

	///
	/// \brief Iterative message passing over the join graph.
	/// \param nIter 	The number of iterations
	/// \param stopTime	The time limit
	/// \param stopObj 	The error tolerance (ie, difference between objective values)
	///
	void propagate(size_t nIter, double stopTime = -1, double stopObj = -1);
};

} // namespace


#endif /* IBM_MERLIN_IJGP_H_ */
