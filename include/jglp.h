/*
 * jglp.h
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

/// \file jglp.h
/// \brief Join Graph Linear Programming (JGLP) algorithm
/// \author Radu Marinescu radu.marinescu@ie.ibm.com


#ifndef IBM_MERLIN_JGLP_H_
#define IBM_MERLIN_JGLP_H_

#include "algorithm.h"
#include "graphical_model.h"

namespace merlin {

/**
 * Join-Graph Linear Programming (ie, Weighted mini-buckets for MAP inference)
 *
 * Tasks supported: MAP
 *
 * JGLP is a specialization of the weighted mini-bucket algorithm to the MAP
 * inference task. It builds a join graph by running the mini-bucket algorithm
 * schemetically, and then it connects the mini-buckets residing in the same
 * bucket. JGLP shifts costs along the edges of the join graph by matching the
 * max marginals of any two clusters connected by an edge. In practice,
 * the algorithm typically tightens the MAP upper bound, however this is not
 * a guarantee. It keeps track of the tightest upper bound found and its
 * corresponding assignment.
 */

class jglp: public graphical_model, public algorithm {
public:
	typedef graphical_model::findex findex;        ///< Factor index
	typedef graphical_model::vindex vindex;        ///< Variable index
	typedef graphical_model::flist flist;          ///< Collection of factor indices
	typedef flist::const_iterator flistIt;		   ///< Iterator for collection of factor indeces

public:

	///
	/// \brief Default constructor.
	///
	jglp() :
		graphical_model() {
		set_properties();
	};

	///
	/// \brief Constructor.
	///
	jglp(const graphical_model& gm) :
		graphical_model(gm), m_gmo(gm) {
		clear_factors();
		set_properties();
	};

	///
	/// \brief Clone the algorithm.
	/// \return the pointer to the object containing the cloned algorithm.
	///
//	virtual jglp* clone() const {
//		jglp* gm = new jglp(*this);
//		return gm;
//	}

	double ub() const {
		return m_logz;
	}
	double lb() const {
		return m_lb;
	}
	std::vector<index> best_config() const {
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
		//return _beliefs[0];
		throw std::runtime_error("Not implemented");
	}
	const factor& belief(variable v) const {
		//return _beliefs[0];
		throw std::runtime_error("Not implemented");
	}
	const factor& belief(variable_set vs) const {
		//return _beliefs[0];
		throw std::runtime_error("Not implemented");
	}
	const std::vector<factor>& beliefs() const {
		//return _beliefs;
		throw std::runtime_error("Not implemented");
	}

	///
	/// \brief Return the original graphical model.
	///
	const graphical_model& get_gm_orig() const {
		return m_gmo;
	}

	///
	/// \brief Properties of the algorithm.
	///
	MER_ENUM( Property , iBound,Order,Iter,Debug );

protected:
	// Members:

	graphical_model m_gmo; 					///< Original graphical model
	OrderMethod m_order_method;				///< Variable ordering heuristic
	size_t m_ibound;						///< Mini-bucket i-bound
	double m_logz;							///< Log partition function
	double m_lb;							///< Lower bound
	variable_order_t m_order; 				///< Elimination ordering
	std::vector<vindex> m_parents; 			///< Pseudo tree
	std::vector<index> m_best_config;		///< MAP assignment
	size_t m_num_iter;						///< Number of iterations
	std::vector<flist> m_mini_buckets;		///< Mini-buucket partitionings
	bool m_debug;							///< Internal debugging flag

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
	void set_pseudo_tree(const std::vector<vindex>& p) {
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
	void set_graphical_model(const std::vector<factor>& fs) {
		m_gmo = graphical_model(fs);
	}

	///
	/// \brief Set the properties of the algorithm.
	/// \param opt 	The string containing comma separated property value pairs
	///
	virtual void set_properties(std::string opt = std::string()) {
		if (opt.length() == 0) {
			set_properties("iBound=4,Order=MinFill,Iter=100,Debug=0");
			return;
		}
		m_debug = false;
		std::vector<std::string> strs = split(opt, ',');
		for (size_t i = 0; i < strs.size(); ++i) {
			std::vector<std::string> asgn = split(strs[i], '=');
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
				m_num_iter = atol(asgn[1].c_str());
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
		return F.max(vs);
	}

	///
	/// \brief Compute the marginal over a set of variables.
	/// \param F 	The reference of the factor to marginalize over
	/// \param vs 	The set of variables representing the scope of the marginal
	/// \return the factor representing the marginal over the set of variables.
	///
	factor marg(const factor& F, const variable_set& vs) {
		return F.maxmarginal(vs);
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
	double score(const std::vector<factor>& fin, const variable& VX, size_t i, size_t j) {
		double err;
		const factor& F1 = fin[i], &F2 = fin[j];           // (useful shorthand)
		size_t iBound = std::max(std::max(m_ibound, F1.nvar() - 1),
				F2.nvar() - 1);      // always OK to keep same size
		variable_set both = F1.vars() + F2.vars();
		if (both.nvar() > iBound+1)
			err = -3;  // too large => -3
		else // greedy scope-based 2 (check if useful???)
			err = 1.0 / (F1.nvar() + F2.nvar());
//		else // scope-based => constant score
//			err = 1;
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
		friend std::ostream& operator<<(std::ostream& out, const sPair& p) {
			out << "(" << p.first << ", " << p.second << ")";
			return out;
		}
		;
	};

	///
	/// \brief Initialize the JGLP algorithm.
	///
	void init();

	///
	/// \brief Join graph propagation for tightening the bound.
	/// \param nIter 	The number of iterations
	/// \param stopTime The time limit
	/// \param stopObj 	The difference between two successive objective values
	///
	void tighten(size_t nIter, double stopTime = -1, double stopObj = -1);

	///
	/// \brief Write the solution to the output file.
	/// \param filename 		The output file name
	/// \param evidence 		The evidence variable value pairs
	/// \param old2new			The mapping between old and new variable indexing
	/// \param orig 			The graphical model prior to asserting evidence
	/// \param output_format	The output format (json or uai)
	///
	void write_solution(std::ostream& out, const std::map<size_t, size_t>& evidence,
			const std::map<size_t, size_t>& old2new, const graphical_model& orig,
			const std::set<size_t>& dummies, int output_format);

	///
	/// \brief Run the JGLP algorithm.
	///
	void run();


	///
	/// \brief Return the MAP configuration following the propagation.
	///
	std::vector<index> config();

private:

	// Partition the bucket into mini-buckets; initially 'ids' is the bucket,
	// but at the end it contains a list of mini-buckets (factors in each
	// mini-bucket are combined into a single one, thus a mini-bucket is a factor
	// defined over at most iBound distinct variables).
	void partition(variable VX, std::vector<factor>& fin, std::vector<flist>& vin,
			std::vector<double>& Norm, std::vector<flist>& Orig,
			std::vector<flist>& New, flist& ids);

};

} // namespace

#endif // re-include
