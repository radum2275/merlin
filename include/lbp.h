/*
 * lbp.h
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

/// \file lbp.h
/// \brief Loopy Belief Propagation (LBP) algorithm
/// \author Radu Marinescu radu.marinescu@ie.ibm.com
///

#ifndef IBM_MERLIN_LOOPY_BP_H_
#define IBM_MERLIN_LOOPY_BP_H_

#include "algorithm.h"
#include "factor_graph.h"
#include "indexed_heap.h"

namespace merlin {

/**
 * Loopy belief propagation over the factor graph.
 *
 * Tasks supported: MAR
 *
 * LBP propagates the messages along a factor graph representation of the
 * original graphical model. It also computes an estimate of the partition
 * function (it doesn't guarantee an upper or lower bound on Z).
 *
 */

class lbp: public algorithm, public factor_graph {
public:
	typedef factor_graph::findex findex;	///< Factor index
	typedef factor_graph::vindex vindex;	///< Variable index
	typedef factor_graph::flist flist;		///< Collection of factor indices

public:
	// Constructors : from nothing, copy, list of factors, or input iterators:

	///
	/// \brief Constructs algorithm instance over an empty factor graph.
	///
	lbp() :	factor_graph() {
		set_properties();
	}

	///
	/// \brief Constructs algorithm instance from a copy of a factor graph.
	/// \param fg	A factor graph
	///
	lbp(const factor_graph& fg) : factor_graph(fg) {
		set_properties();
	}

	///
	/// \brief Constructs algorithm instance from a list of factors.
	/// \param fs 	A list of factors
	///
	lbp(std::vector<factor> fs) : factor_graph(fs) {
		set_properties();
	}

	///
	/// \brief Constructs algorithm instance from input iterators.
	/// \param f 	An iterator to beginning
	/// \param l	An iterator to end
	///
	template<class InputIterator>
	lbp(InputIterator f, InputIterator l) :	factor_graph(f, l) {
		set_properties();
	}

	///
	/// \brief Clone the algorithm object.
	/// \return the pointer to the object containing the cloned algorithm.
	///
//	virtual lbp* clone() const {
//		lbp* fg = new lbp(*this);
//		return fg;
//	}

	///
	/// \brief Mutable accessor to a belief.
	///
	merlin::factor& bel(size_t f) {
		return m_beliefs[f];
	}

	const factor& belief(size_t f) const {
		return m_beliefs[f];
	}
	const factor& belief(variable v) const {
		return belief(local_factor(v));
	}
	const factor& belief(variable_set vs) const {
		throw std::runtime_error("Not implemented");
	}
	const std::vector<factor>& beliefs() const {
		return m_beliefs;
	}

	// Not a bound-producing algorithm but can try to produce a good config:

	double lb() const {
		throw std::runtime_error("Not available");
	}
	double ub() const {
		throw std::runtime_error("Not available");
	}
	std::vector<index> best_config() const {
		throw std::runtime_error("Not available");
	}

	// Gives an estimate of the partition function, but not a bound:

	double logZ() const {
		return m_logz;
	}
	double logZub() const {
		throw std::runtime_error("Not available");
	}
	double logZlb() const {
		throw std::runtime_error("Not available");
	}

	///
	/// \brief Types of propagation schedules.
	///
	MER_ENUM( Schedule , Fixed,Random,Flood,Priority );

	///
	/// \brief Properties of the algorithm.
	///
	MER_ENUM( Property , Schedule,Distance,StopIter,StopObj,StopMsg,Debug );

	///
	/// \brief Set the properties of the algorithm.
	/// \param opt 	The string containing comma separated property value pairs
	///	
	virtual void set_properties(std::string opt = std::string()) {
		if (opt.length() == 0) {
			set_properties(
					"Schedule=Priority,Distance=HPM,StopIter=10,StopObj=-1,StopMsg=-1,Debug=0");
			return;
		}
		m_debug = false;
		std::vector<std::string> strs = split(opt, ',');
		for (size_t i = 0; i < strs.size(); ++i) {
			std::vector<std::string> asgn = split(strs[i], '=');
			switch (Property(asgn[0].c_str())) {
			case Property::Schedule:
				m_sched = Schedule(asgn[1].c_str());
				break;
			case Property::Distance:
				m_dist = factor::Distance(asgn[1].c_str());
				break;
			case Property::StopIter:
				set_stop_iter(strtod(asgn[1].c_str(), NULL));
				break;
			case Property::StopObj:
				set_stop_obj(strtod(asgn[1].c_str(), NULL));
				break;
			case Property::StopMsg:
				set_stop_msg(strtod(asgn[1].c_str(), NULL));
				break;
			case Property::Debug:
				m_debug = (atol(asgn[1].c_str()) == 0) ? false : true;
				break;
			default:
				break;
			}
		}
	}

	// Initialize the data structures:

	///
	/// \brief Initialize the LBP algorithm.
	///
	void init();

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


	// Run algorithm LBP
	void run();

protected:

	// Contained objects
	std::vector<factor> m_beliefs;      ///< Store calculated messages and beliefs
	std::vector<factor> m_msg;			///< Store messages
	std::vector<factor> m_msg_new;		///< Store new messages
	indexed_heap m_priority;            ///< Store m_priority schedule of edges
	std::vector<findex> m_forder;	    ///< Fixed order of factors
	Schedule m_sched;                   ///< Schedule type
	factor::Distance m_dist;	        ///< Message distance measure for priority_
	double m_logz;                      ///< Current objective function value
	bool m_debug;						///< Internal debugging flag

	/// 
	/// \brief Entropy computation.
	///
	/// Calculate the entropy contribution to the free energy from a node
	/// \param n 	The index of the node
	/// \return the entropy contribution from node *n*.
	double obj_entropy(size_t n) {
		double obj = belief(n).entropy();
		if (!is_var_node(n)) {
			variable_set vs = adjacent_vars(n);
			for (variable_set::const_iterator i = vs.begin(); i != vs.end(); ++i)
				obj -= belief(n).marginal(*i).entropy();
		}
		return obj;
	}

	///
	/// \brief Belief computation.
	///
	/// Re-calculate the belief at a node from the current incoming messages.
	/// \param n 	The index of the node
	///
	void calc_belief(size_t n) {
		const set<edge_id>& nbrs = neighbors(n);		// get all incoming edges
		bel(n) = factor(n);             // calculate local factor times messages
		for (set<edge_id>::const_iterator i = nbrs.begin(); i != nbrs.end(); ++i)
			bel(n) *= m_msg[i->ridx];
		bel(n) /= bel(n).sum();                 // and normalize
	}

	///
	/// \brief Incoming messages computation.
	///
	/// Accept all the incoming messages into a node, and recompute its belief.
	/// \param n 	The index of the node
	///
	void accept_incoming(size_t n) {
		const set<edge_id>& nbrs = neighbors(n); // get the list of neighbors
		bel(n) = get_factor(n); //   and start with just the local factor
		for (set<edge_id>::const_iterator i = nbrs.begin();
				i != nbrs.end(); ++i) {
			m_msg[i->ridx] = m_msg_new[i->ridx]; // accept each new incoming message
			bel(n) *= m_msg[i->ridx];           //   and include it in the belief
			if (m_sched == Schedule::Priority) {
				m_priority.erase(i->ridx); // accepted => remove from priority_ queue
			}
		}
		bel(n) /= bel(n).sum(); // normalize belief
	}

	///
	/// \brief Outgoing (new) messages computation.
	///
	/// Recompute new messages from node n to its neighbors.
	/// \param n 	The index of the node
	///
	void update_outgoing(size_t n) {
		const set<edge_id>& nbrs = neighbors(n); // get the list of neighbors
		for (set<edge_id>::const_iterator i = nbrs.begin();
				i != nbrs.end(); ++i) {
			m_msg_new[i->idx] = (belief(n) / m_msg[i->ridx]).marginal(belief(i->second).vars());
			m_msg_new[i->idx] /= m_msg_new[i->idx].sum();   // normalize message
			if (m_sched == Schedule::Priority) { // and update m_priority in schedule
				m_priority.insert(m_msg_new[i->idx].distance(m_msg[i->idx], m_dist),
						i->idx);
			}
		}
	}

};

} // namespace

#endif /* IBM_MERLIN_LBP_H_ */
