/*
 * cte.h
 *
 *  Created on: 14 Sep 2018
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

/// \file cte.h
/// \brief Clique-Tree Elimination for joint marginals
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_MERLIN_CTE_H_
#define IBM_MERLIN_CTE_H_


#include "graphical_model.h"
#include "algorithm.h"
#include "observation.h"

namespace merlin {

/**
 * Clique-Tree Elimination (CTE)
 *
 * Tasks supported: MAR and joint MAR
 *
 * Clique-Tree Elimination is an exact inference algorithm that supports the
 * main inference tasks in graphical models (PR, MAR, MAP, MMAP). It belongs to
 * the class on join-tree elimination algorithms and is driven by a clique-tree
 * structure. Its complexity is time and space exponential in the induced width.
 * Note that this is meant to be a reference implementation and therefore is
 * applicable to relatively small models having relatively small induced widths.
 */

namespace detail {
	class edge;

	class node {
	public:
		size_t id;
		variable_set clique;
		factor theta;
		node* parent;
		vector<node*> children;
		double weight;
		vector<graphical_model::findex> originals;
		factor belief;
		vector<edge*> edges;
		node() {
			clear();
		};
		void clear() {
			clique.clear();
			theta = belief = factor(1.0);
			parent = NULL; children.clear();
			weight = 0.0;
		}
	};

	class edge {
	public:
		variable_set sepset;
		node* first;
		node* second;
		factor fwd; // message from 1 to 2
		factor bwd; // message from 2 to 1
		void messageFwd();
		void messageFwd(std::vector<int>& evidence);
		void messageBwd();
		void messageBwd(std::vector<int>& evidence);
		edge() {
			clear();
		}
		edge(node* ni, node* nj) {
			first = ni;
			second = nj;
			sepset = (ni->clique & nj->clique);
			fwd = bwd = factor(sepset, 1.0);
		}
		void clear() {
			sepset.clear();
			fwd = bwd = factor(1.0);
			first = second = NULL;
		}
		void reset() {
			fwd = bwd = factor(sepset, 1.0);
		}
	};

} // namespace


class cte: public graphical_model, public algorithm {
public:
	typedef graphical_model::findex findex;        ///< Factor index
	typedef graphical_model::vindex vindex;        ///< Variable index
	typedef graphical_model::flist flist;          ///< Collection of factor indices

public:

	///
	/// \brief Default constructor.
	///
	cte() : graphical_model() {
		set_properties();
	}

	///
	/// \brief Constructor with a graphical model.
	///
	cte(const graphical_model& gm) : graphical_model(gm), m_gmo(gm) {
		clear_factors();
		set_properties();
	}

	~cte() {
		for (vector<detail::edge*>::iterator i = m_edges.begin();
				i != m_edges.end(); ++i) {
			delete *i;
		}
		m_edges.clear();
		m_messages.clear();
	};

	///
	/// \brief Clone the algorithm.
	/// \return the pointer to the new object containing the cloned algorithm.
	///
//	inline virtual cte* clone() const {
//		cte* gm = new cte(*this);
//		return gm;
//	}

	///
	inline double ub() const {
		return m_logz;
	}
	inline double lb() const {
		throw std::runtime_error("Not implemented");
	}
	inline std::vector<size_t> best_config() const {
		throw std::runtime_error("Not implemented");
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

	inline const factor& belief(size_t f) const {
		return m_beliefs[f];
	}
	inline const factor& belief(variable v) const {
		return m_beliefs[v];
	}
	inline const factor& belief(variable_set vs) const {
		throw std::runtime_error("Joint belief must be contained into exactly one cluster.");
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
	/// \brief Run the clique-tree elimination algorithm.
	///
	void run();

	///
	/// \brief Inference tasks supported.
	///
	MER_ENUM( Task, PR,MAR,MAP,MMAP );

	///
	/// \brief Properties of the algorithm
	///
	MER_ENUM( Property , Order,Task,Debug,Verbose );


	// Setting properties (directly or through property string):

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
			set_properties("Order=MinFill,Task=MAR,Debug=0,Verbose=1");
			return;
		}
		m_debug = false;
		m_verbose = 1;
		std::vector<std::string> strs = merlin::split(opt, ',');
		for (size_t i = 0; i < strs.size(); ++i) {
			std::vector<std::string> asgn = merlin::split(strs[i], '=');
			switch (Property(asgn[0].c_str())) {
			case Property::Order:
				m_order.clear();
				m_order_method = graphical_model::OrderMethod(asgn[1].c_str());
				break;
			case Property::Task:
				m_task = Task(asgn[1].c_str());
				break;
			case Property::Debug:
				m_debug = (atol(asgn[1].c_str()) == 0) ? false : true;
				break;
			case Property::Verbose:
				m_verbose = atol(asgn[1].c_str());
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
	/// \brief Eliminate a set of variables from a factor.
	/// \param F 	The reference of the factor to eliminate from
	///	\param vs 	The set of variables to be eliminated
	/// \param w 	The weight of the weighted elimination operator
	/// \return the factor resulted from eliminating the set of variables.
	///
	inline factor elim(const factor& F, const variable_set& vs) {
		return F.sum(vs);
	}

	///
	/// \brief Compute the weighted marginal over a set of variables.
	/// \param F 	The reference of the factor to marginalize over
	/// \param vs 	The set of variables representing the scope of the marginal
	/// \param w 	The weight of the weighted elimination operator
	/// \return the factor representing the weighted marginal over the set of variables.
	///
	inline factor marg(const factor& F, const variable_set& vs) {
		return F.marginal(vs);
	}

	///
	/// \brief Initialize the clique tree elimination algorithm.
	///
	void init();

	///
	/// \brief Reinitialize the clique tree with new original factors
	///
	void reinit(const std::vector<factor>& factors);

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
	/// \brief Forward (top-down) message passing along the edge of the bucket tree.
	///
	void forward();
	void forward(std::vector<int>& evidence);

	///
	/// \brief Backward (bottom-up) message passing along the edges of the bucket tree.
	///
	void backward();
	void backward(std::vector<int>& evidence);

	///
	/// \brief Propagate the messages along the edges of the bucket tree.
	///
	void calibrate();

	///
	/// \brief Update the beliefs (marginals or max-marginals)
	///
	void update();

	///
	/// \brief Propagate the evidence in the join tree
	/// \param evidence	The evidence vector
	///
	bool propagate_evidence(std::vector<int>& evidence);

	///
	/// \brief Compute the joint marginal of a set of variables
	/// \param scope	The variables of the joint marginal
	///
	void joint_marginal(const variable_set& scope);

	///
	/// \brief Compute the joint marginal of a set of variables subject to evidence
	/// \param scope	The variables of the joint marginal
	/// \param evidence	The evidence vector
	///
	void joint_marginal(const variable_set& scope, std::vector<int>& evidence);

	///
	/// \brief Return the joint marginal
	///
	inline const factor& get_joint_marginal() {
		return m_marginal;
	}

private:

	///
	/// \brief Build the clique tree from a set of maximal cliques
	///
	void build_clique_tree(std::vector<std::set<size_t> >& cliques);

	///
	/// \brief Check if a configuration is compatible with an evidence vector
	/// \param config	The configuration as a map of variable value pairs
	/// \param evidence	The evidence vector
	///	\return true if all variables in the configuration have the same values
	///				 as in the evidence vector
	///
	bool is_compatible(std::map<size_t, size_t>& config,
			std::vector<int>& evidence);

protected:
	// Members:

	graphical_model m_gmo; 				///< Original graphical model
	Task m_task;						///< Inference task
	OrderMethod m_order_method;			///< Variable ordering method
	double m_logz;						///< Log partition function value
	variable_order_t m_order;			///< Variable order
	vector<factor> m_beliefs; 			///< Marginals
	std::vector<vindex> m_query; 		///< Variables for the joint marginal
	factor m_marginal;					///< Joint marginal of the query variables
	std::vector<int> m_evidence;		///< Evidence propagated during EM learning

private:
	// CT local structures:

	detail::node* m_root;				///< Clique tree root
	vector<detail::node> m_clusters;	///< Cliques
	vector<detail::edge*> m_edges;		///< Edges
	vector<detail::edge*> m_messages;	///< Propagation schedule
	vector<int> m_var2clique;			///< Variable to clique map

	bool m_debug;						///< Internal flag for debugging only
	int m_verbose;						///< Verbosity level (default = 1)
};

} // namespace




#endif /* IBM_MERLIN_CTE_H_ */
