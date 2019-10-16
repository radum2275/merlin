/*
 * factor_graph.h
 *
 *  Created on: Feb 8, 2013
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

/// \file factor_graph.h
/// \brief A factor graph graphical model
/// \author Radu Marinescu radu.marinescu@ie.ibm.com


#ifndef IBM_MERLIN_FACTORGRAPH_H_
#define IBM_MERLIN_FACTORGRAPH_H_

#include "graphical_model.h"

namespace merlin {

///
/// \brief A factor graph base class
///
/// A graphical model represented as a bipartite graph between *variable nodes*,
/// corresponding to the variables, and *factor nodes*, corresponding to the factors.
/// Internally, factors and variables are mainly referenced by integer indices, with
/// 0 <= f < nFactors() and 0 <= v < nvar(). However, many interface functions are called
/// using a variable object, var(label,dim). To convert, var(v) gives the vth var object
/// and the internal function _vindex(V) gives the index corresponding to variable object V.
///
class factor_graph: public graphical_model {
public:
	typedef graphical_model::findex findex;		///< Factor index
	typedef graphical_model::vindex vindex;		///< Variable index
	typedef graphical_model::flist flist; 		///< Collection of factor indices

protected:
	typedef size_t eindex;						///< Edge index

public:

	///
	/// \brief Default constructor.
	///
	factor_graph() :
		graphical_model(), m_vindex() {
	};

	///
	/// \brief Copy constructor.
	/// \param fg 	The factor graph to be copied
	///
	factor_graph(const factor_graph& fg) :
			graphical_model((graphical_model&) fg), m_vindex(fg.m_vindex) {
	}

	///
	/// \brief Constructor.
	/// \param fs 	The list of factors
	///
	factor_graph(std::vector<factor> fs) : graphical_model(fs), m_vindex() {
		m_vindex.resize(nvar());
		create_factor_graph();
	}

	///
	/// \brief Constructor.
	/// \param first 	The *begin* iterator of the factors list
	/// \param last 	The *end* iterator of the factors list
	///
	template<class InputIterator>
	factor_graph(InputIterator first, InputIterator last) :
			graphical_model(first, last), m_vindex() {
		m_vindex.resize(nvar());
		create_factor_graph();
	}

	///
	/// \brief Clone the factor graph.
	/// \return a pointer to the cloned factor graph.
	///
	virtual factor_graph* clone() {
		factor_graph* fg = new factor_graph(*this);
		return fg;
	}


	///
	/// \brief Add a new factor to the graphical model.
	/// \param F 	The factor to be added
	/// \return the index of the newly added factor.
	findex add_factor(const factor& F) {
		findex use = graphical_model::add_factor(F); // add the factor to the underlying collection
		m_vindex.resize(nvar(), vindex(-1));	// then update variable nodes and edges (!!!)
		if (F.nvar() == 1 && m_vindex[*F.vars().begin()] == vindex(-1))
			m_vindex[*F.vars().begin()] = use;
		else
			for (variable_set::const_iterator v = F.vars().begin();
					v != F.vars().end(); ++v) {
				if (m_vindex[*v] == vindex(-1))
					m_vindex[*v] = graphical_model::add_factor(factor(*v, 1.0));
				add_edge(use, local_factor(*v));		// add edge to variable node
			}
		return use;
	}

	///
	/// \brief Remove a factor from the model.
	/// \param f 	The index of the factor to be removed
	///
	void remove_factor(findex f) {
		variable_set vs = m_factors[f].vars();		// save the variables for clean-up
		graphical_model::remove_factor(f);// remove the factor from the collection
		for (variable_set::const_iterator v = vs.begin(); v != vs.end(); ++v) {
			remove_edge(f, local_factor(*v));	// remove the edge to the variable node
			if (local_factor(*v) == f)
				m_vindex[_vindex(*v)] = vindex(-1);// removal of variable node => mark as missing (!!!)
		}
	}

	///
	/// \brief Retrieve the factor index corresponding to a variable node.
	/// \param i 	The index of a variable
	/// \return the index of the factor associated with the variable in
	/// 	the factor graph.
	///
	findex local_factor(vindex i) const {
		return m_vindex[i];
	}

	///
	/// \brief Retrieve the factor index corresponding to a variable node.
	/// \param v 	The index of a variable
	/// \return the index of the factor associated with the variable in
	/// 	the factor graph.
	///	
	findex local_factor(variable v) const {
		return m_vindex[_vindex(v)];
	}

	///
	/// \brief Check if a factor is a variable node in the factor graph.
	/// \param i 	The index of the factor
	/// \return *true* if the factor corresponds to a variable node, 
	/// 	and *false* otherwise.
 	///
	bool is_var_node(findex i) const {
		return (get_factor(i).nvar() == 1
				&& local_factor(*get_factor(i).vars().begin()) == i);
	}

	///
	/// \brief Retrieve the factors adjacent to a variable node in the factor graph.
	/// \param v 	The index of the variable
	/// \return the list of factor indexes corresponding to the adjacent factors.
	///
	flist adjacent_factors(variable v) const {
		return _neighbors(local_factor(v));
	}

	///
	/// \brief Retrieve the factors adjacent to a variable node in the factor graph.
	/// \param v 	The index of the variable
	/// \return the list of factor indexes corresponding to the adjacent factors.
	///	
	flist adjacent_factors(vindex v) const {
		return _neighbors(local_factor(v));
	}

	///
	/// \brief Retrieve the variables adjacent to a factor in the factor graph.
	/// \param f 	The index of the factor
	/// \return the list of variable indexes corresponding to the adjacent variables.
	///	
	variable_set adjacent_vars(findex f) const {
		return get_factor(f).vars();
	}

	///
	/// \brief Swap the contents of two factor graphs.
	/// \param gm 	The factor graph to swap with
	///
	void swap(factor_graph& gm);

protected:

	// Functions:

	///
	/// \brief Retrieve the neighbors of a factor.
	/// \param i 	The index of the factor
	/// \return the list of factor indexes that are connected to the input factor.
	///
	flist _neighbors(findex i) const {
		flist fl;
		for (set<edge_id>::const_iterator it = neighbors(i).begin();
				it != neighbors(i).end(); ++it)
			fl += it->second;
		return fl;
	}


	///
	/// \brief Create the factor graph.
	///
	/// Creates the bi-partite graph by adding nodes corresponding to local
	/// variable factors (ie, unary factors definded for each variable)
	/// and connecting them to the nodes corresponding to the factors. Let M be
	/// the number of factors. Then the nodes corresponding to the factors
	/// are indexed [0 .. M-1], while the nodes corresponding to the variables
	/// are indexed [M .. M+N-1] where N is the number of variables.
	///
	void create_factor_graph() {

		// init var to node index map
		if (m_vindex.size() < nvar()) {
			m_vindex.resize(nvar());
		}

		// identify unary factors (if any)
		std::vector<bool> found(nvar(), false);
		for (size_t i = 0; i < m_factors.size(); ++i) {
			if (m_factors[i].nvar() == 1) {
				size_t v = _vindex(*m_factors[i].vars().begin());
				if (!found[v]) m_vindex[v] = i;
				found[v] = true;
			}
		}

		// create local factors for each of the variables
		for (size_t v = 0; v < found.size(); ++v) {
			if (!found[v]) {
				 m_vindex[v] = graphical_model::add_factor(factor(var(v), 1.0));
			}
		}

		// create the bi-partite graph (edges from factor node to its variable nodes)
		if (num_edges() != 0) {
			throw std::runtime_error("Initial factor graph must be empty.");
		}
		for (size_t i = 0; i < m_factors.size(); ++i) {
			if (!is_var_node(i)) { // factor node
				for (variable_set::const_iterator v = m_factors[i].vars().begin();
						v != m_factors[i].vars().end(); ++v) {
					add_edge(i, local_factor(*v));
				}
			}
		}
	}
	
	///
	/// \brief Retrieve the index of the edge between two nodes.
 	/// \param i 	The (factor) index of the head
 	/// \param j 	The (factor) index of the tail
 	/// \return the index of the edge between the two factors.
	eindex _eindex(findex i, findex j) {
		return edge(i, j).idx;
	}

	// Queue type wrappers for abstracting depth-1st, width-1st, etc. 
	// for tree search:

	///
	/// \brief Queue type wrapper.
	///
	template<class T> class abstract_queue {
	public:
		virtual ~abstract_queue() {
		}
		;
		virtual void pop()=0;
		virtual void push(T)=0;
		virtual T first()=0;
		virtual bool empty()=0;
	};

	///
	/// \brief Stack type for depth-first search.
	///
	template<class T> class abstract_queue_lifo: public abstract_queue<T>,
			public std::stack<T> {
	public:
		void push(T t) {
			std::stack<T>::push(t);
		}
		void pop() {
			std::stack<T>::pop();
		}
		T first() {
			return this->top();
		}
		bool empty() {
			return std::stack<T>::empty();
		}
	};
	//using std::stack<T>::push; using std::stack<T>::pop; T first() { return this->top(); } };

	///
	/// \brief Queue type for breadth-first search.
	///
	template<class T> class abstract_queue_fifo: public abstract_queue<T>,
			public std::queue<T> {
	public:
		void push(T t) {
			std::queue<T>::push(t);
		}
		void pop() {
			std::queue<T>::pop();
		}
		T first() {
			return this->front();
		}
		bool empty() {
			return std::queue<T>::empty();
		}
	};

public:

	///
	/// \brief Type of tree.
	///
	MER_ENUM( TreeType , WidthFirst,DepthFirst,MaxWeight );

	///
	/// \brief Create a spanning tree of the factor graph.
	/// \param tt 	The type of spanning tree
	/// \param root The root (variable index) of the spanning tree
	/// \return the spanning tree as a list of edge indexes in the factor graph.
	///
	std::vector<edge_t> span_tree(TreeType tt, variable root) {
		abstract_queue<findex> *Q;
		switch (tt) {
		case TreeType::WidthFirst:
			Q = new abstract_queue_fifo<findex>();
			break;
		case TreeType::DepthFirst:
			Q = new abstract_queue_lifo<findex>();
			break;
		case TreeType::MaxWeight:
		default:
			throw std::runtime_error("Not implemented");
			break;
		}
		Q->push(local_factor(root));
		std::vector<findex> used, fUsed;
		used.resize(nvar(), findex(-1));
		fUsed.resize(num_factors(), 0);
		std::vector<edge_t> tree;
		tree.reserve(2 * nvar());
		while (!Q->empty()) {
			findex next = Q->first();
			Q->pop();
			if (fUsed[next] != 0)
				continue;						// already taken this factor?
			size_t nFound = 0, parent = 0;// if not look over variables and make sure
			const variable_set& vs = get_factor(next).vars();	//   at most one is already covered;
			for (variable_set::const_iterator v = vs.begin(); v != vs.end(); ++v) {
				if (used[_vindex(*v)] != vindex(-1)) {
					parent = used[_vindex(*v)];
					nFound++;
				}
			}						// if one is covered that's our parent; more
			if (nFound > 1)
				continue;//   means this factor is unavailable so keep looking
			fUsed[next] = 1;	// otherwise, mark factor, variables as covered
			if (is_var_node(next))
				used[_vindex(*vs.begin())] = next;// always cover by var node once visited
			else
				for (variable_set::const_iterator v = vs.begin(); v != vs.end();
						++v) {
					if (used[_vindex(*v)] == vindex(-1))
						used[_vindex(*v)] = next;// otherwise remember which factor we are
				}
			// if we have a parent, add an edge to our trees
			if (nFound > 0)
				tree.push_back(edge_t(parent, next)); //, _eindex(parent,next),_eindex(next,parent))); //!!!
			// if this is a root, do nothing (!!!); don't keep track of isolated nodes?

			const set<edge_id> nbrs = neighbors(next);// add all neighbors to the search queue
			for (set<edge_id>::const_iterator n = nbrs.begin(); n != nbrs.end();
					++n) {
				Q->push(n->second);
			}
		}
		delete Q;
		return tree;
	}

	///
	/// \brief Create a spanning tree of the factor graph rooted by a randomly selected variable.
	/// \param tt 	The type of spanning tree
	/// \return the spanning tree as a list of edge indexes in the factor graph.
	///
	std::vector<edge_t> span_tree(TreeType tt = TreeType::WidthFirst) {
		return span_tree(tt, var(randi(nvar())));
	}

protected:
	// Contained objects

	std::vector<size_t> m_vindex;	///< Factors representing variable nodes in factor graph

};

} // namespace

#endif  // re-include
