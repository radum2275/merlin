/*
 * lbp.h
 *
 *  Created on: 24 Mar 2015
 *      Author: radu
 */

/// \file lbp.h
/// \brief Loopy Belief Propagation (LBP) algorithm
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
	lbp() :
			factor_graph() {
		set_properties();
	}

	///
	/// \brief Constructs algorithm instance from a copy of a factor graph.
	/// \param fg	A factor graph
	///
	lbp(const factor_graph& fg) :
			factor_graph(fg) {
		set_properties();
	}

	///
	/// \brief Constructs algorithm instance from a list of factors.
	/// \param fs 	A list of factors
	///
	lbp(std::vector<factor> fs) :
			factor_graph(fs) {
		set_properties();
	}

	///
	/// \brief Constructs algorithm instance from input iterators.
	/// \param f 	An iterator to beginning
	/// \param l	An iterator to end
	///
	template<class InputIterator>
	lbp(InputIterator f, InputIterator l) :
			factor_graph(f, l) {
		set_properties();
	}

	///
	/// \brief Clone the algorithm object.
	/// \return the pointer to the object containing the cloned algorithm.
	///
	virtual lbp* clone() const {
		lbp* fg = new lbp(*this);
		return fg;
	}

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
		return m_log_z;
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
	MER_ENUM( Property , Schedule,Distance,stopIter,stopObj,stopMsg );

	///
	/// \brief Set the properties of the algorithm.
	/// \param opt 	The string containing comma separated property value pairs
	///	
	virtual void set_properties(std::string opt = std::string()) {
		if (opt.length() == 0) {
			set_properties(
					"Schedule=Priority,Distance=HPM,stopIter=10,stopObj=-1,stopMsg=-1");
			return;
		}
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
			case Property::stopIter:
				set_stop_iter(strtod(asgn[1].c_str(), NULL));
				break;
			case Property::stopObj:
				set_stop_obj(strtod(asgn[1].c_str(), NULL));
				break;
			case Property::stopMsg:
				set_stop_msg(strtod(asgn[1].c_str(), NULL));
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
	void init() {
		// Start the timer and store it
		m_start_time = timeSystem();

		// Prologue
		std::cout << VERSIONINFO << std::endl << COPYRIGHT << std::endl;
		std::cout << "Initialize inference engine ..." << std::endl;
		std::cout << "+ tasks supported  : PR,MAR" << std::endl;
		std::cout << "+ algorithm        : " << "LBP" << std::endl;
		std::cout << "+ inference task   : " << "MAR" << std::endl;
		std::cout << "+ schedule         : " << m_sched << std::endl;

		m_beliefs = std::vector<factor>(m_factors); // copy initial beliefs from factors
		m_msg = std::vector<factor>();
		m_msg.resize(2 * num_edges());        // initialize messages to the identity
		for (size_t e = 0; e < 2 * num_edges(); ++e)
			if (edge(e) != edge_id::NO_EDGE) {  // f'n of the right variables
				m_msg[e] = merlin::factor(
						merlin::factor(edge(e).first).vars()
								& factor(edge(e).second).vars(), 1.0);
			}
		m_msg_new = std::vector<merlin::factor>(m_msg); // copy that as "updated" message list

		m_log_z = 0.0;                    // compute initial partition f'n estimate
		for (size_t f = 0; f < num_factors(); ++f) {
			bel(f) /= bel(f).sum();                     // normalize the beliefs
			m_log_z += (bel(f) * log(factor(f))).sum() + obj_entropy(f); // and compute the free energy estimate
		}

		if (m_sched == Schedule::Priority) {           // for priority scheduling
			for (size_t e = 0; e < 2 * num_edges(); ++e) //  initialize all edges to infinity
				if (edge(e) != edge_id::NO_EDGE)
					m_priority.insert(std::numeric_limits<double>::infinity(), e);
		} else {
			for (size_t f = 0; f < num_factors(); ++f)
				m_forder.push_back(f); // for fixed scheduling, get a default order
		}
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
			throw std::runtime_error("Error while opening output file.");
		}

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
	}

	// Run algorithm LBP
	void run() {
		init();

		size_t stopIter = m_stop_iter * num_factors();// it's easier to count updates than "iterations"

		double dObj = m_stop_obj + 1.0, dMsg = m_stop_msg + 1.0;// initialize termination values
		size_t iter = 0, print = 1; // count updates and "iterations" for printing
		size_t f, n = 0;                                                 //

		std::cout << "Begin message passing over factor graph ..." << std::endl;
		std::cout << "- stopObj  : " << m_stop_obj << std::endl;
		std::cout << "- stopMsg  : " << m_stop_msg << std::endl;
		std::cout << "- stopIter : " << m_stop_iter << std::endl;
		for (; dMsg >= m_stop_msg && iter < stopIter && dObj >= m_stop_obj;) {

			if (m_sched == Schedule::Priority) {          // priority schedule =>
				f = edge(m_priority.top().second).second;	//   get next factor for update from queue
				m_priority.pop();
			} else {                                        // fixed schedule =>
				f = m_forder[n];                    //   get next factor from list
				if (++n == m_forder.size())
					n = 0;
			}

			if (m_sched != Schedule::Flood) {       // For non-"flood" schedules,
				merlin::factor logF = log(factor(f));
				dObj = 0.0;          // compute new belief and update objective:
				dObj -= (belief(f) * logF).sum() + obj_entropy(f); //   remove old contribution
				accept_incoming(f);        //   accept all messages into factor f
				m_log_z += dObj += (belief(f) * logF).sum() + obj_entropy(f); //   re-add new contribution
			}
			update_outgoing(f);		//   update outgoing messages from factor f

			if (m_sched == Schedule::Priority)
				dMsg = m_priority.top().first; // priority schedule => easy to check msg stop
			else if (m_stop_msg > 0 && n == 0) { // else check once each time through all factors
				dMsg = 0.0;
				for (size_t e = 0; e < 2 * num_edges(); ++e)
					dMsg = std::max(dMsg, m_msg_new[e].distance(m_msg[e], m_dist));
			}

			if (m_sched == Schedule::Flood && n == 0) { // for flooding schedules, recalculate all
				dObj = m_log_z;
				m_log_z = 0.0;                   //   the beliefs and objective now
				for (size_t f = 0; f < num_factors(); ++f) {
					accept_incoming(f);
					m_log_z += (belief(f) * log(factor(f))).sum() + obj_entropy(f);
				}
				dObj -= m_log_z;
			}

			if (iter > print * num_factors()) {
				print++;
	//				std::cout << "  LBP: " << _lnZ << " \t\t(dObj=" << dObj << " \t dMsg=" << dMsg << ")\n";

				std::cout << "  LBP: " << std::fixed << std::setw(12) << std::setprecision(6)
					<< m_log_z << " (" << std::scientific << std::setprecision(6)
					<< std::exp(m_log_z) << ") ";
				std::cout << "\td=" << std::scientific << std::setprecision(6)
					<< dObj << "\tm=" << dMsg << "\t time="  << std::fixed
					<< std::setprecision(6)	<< (timeSystem() - m_start_time)
					<< "\ti=" << iter << std::endl;
			}

			iter++;
		}

		// Output solution (UAI output format)
		std::cout << "Converged after " << iter << " iterations in "
			<< (timeSystem() - m_start_time) << std::endl;
		std::cout << "PR" << std::endl;
		std::cout << std::fixed << std::setprecision(6)
			<< m_log_z << " (" << std::scientific << std::setprecision(6)
			<< std::exp(m_log_z) << ")" << std::endl;
		std::cout << "MAR" << std::endl;
		std::cout << nvar();
		for (size_t v = 0; v < m_vindex.size(); ++v) {
			variable VX = var(v);
			std::cout << " " << VX.states();
			for (size_t j = 0; j < VX.states(); ++j) {
				std::cout << " " << std::fixed << std::setprecision(6) << belief(VX)[j];
			}
		}
		std::cout << std::endl;

		std::cout << "Beliefs\n";
		for (size_t i = 0; i < m_beliefs.size(); ++i) {
			std::cout << m_beliefs[i] << std::endl;
		}
	}

protected:

	// Contained objects
	std::vector<factor> m_beliefs;      ///< Store calculated messages and beliefs
	std::vector<factor> m_msg;			///< Store messages
	std::vector<factor> m_msg_new;		///< Store new messages
	indexed_heap m_priority;            ///< Store m_priority schedule of edges
	std::vector<findex> m_forder;	    ///< Fixed order of factors
	Schedule m_sched;                   ///< Schedule type
	factor::Distance m_dist;	        ///< Message distance measure for priority_
	double m_log_z;                     ///< Current objective function value

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
	void accept_incoming(size_t n) {								//
		const set<edge_id>& nbrs = neighbors(n);		// get the list of neighbors
		bel(n) = factor(n);            //   and start with just the local factor
		for (set<edge_id>::const_iterator i = nbrs.begin(); i != nbrs.end();
				++i) {
			m_msg[i->ridx] = m_msg_new[i->ridx]; // accept each new incoming message
			bel(n) *= m_msg[i->ridx];           //   and include it in the belief
			if (m_sched == Schedule::Priority)
				m_priority.erase(i->ridx); // accepted => remove from priority_ queue
		}
		bel(n) /= bel(n).sum();                 // normalize belief
	}

	///
	/// \brief New messages computation.
	///
	/// Recompute new messages from node n to its neighbors.
	/// \param n 	The index of the node
	///
	void update_outgoing(size_t n) {									//
		const set<edge_id>& nbrs = neighbors(n);		// get the list of neighbors
		for (set<edge_id>::const_iterator i = nbrs.begin(); i != nbrs.end();
				++i) {
			m_msg_new[i->idx] = (belief(n) / m_msg[i->ridx]).marginal(
					belief(i->second).vars());
			m_msg_new[i->idx] /= m_msg_new[i->idx].sum();   // normalize message
			if (m_sched == Schedule::Priority) // and update m_priority in schedule
				m_priority.insert(m_msg_new[i->idx].distance(m_msg[i->idx], m_dist),
						i->idx);
		}
	}

};

} // namespace

#endif /* LIB_LBP_H_ */
