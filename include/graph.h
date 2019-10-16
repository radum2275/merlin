/*
 * graph.h
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

/// \file graph.h
/// \brief An undirected graph structure
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_MERLIN_GRAPH_H_
#define IBM_MERLIN_GRAPH_H_

#include "base.h"
#include "set.h"
#include "variable_set.h"

namespace merlin {

typedef std::pair<size_t, size_t> edge_t; ///< Basic Edge type.

///
/// \brief Edge in the graph.
///
/// EdgeID represents the edge (i->j) with index *idx*. The reversed index *ridx*
/// represents the edge (j->i). Blank (or missing) edges can sometimes appear, and
/// they are represented by EdgeID::NO_EDGE (for example, if edge (i,j) is 
/// requested and does not exist).
///
struct edge_id {
	typedef size_t index; 	///< Basic indexing type for the edge.
	index first;			///< Head of the edge.
	index second;			///< Tail of the edge.
	index idx;				///< Index of the edge.
	index ridx;				///< Index of the reversed edge.

	///
	/// \brief Default constructor.
	///
	edge_id() : first(), second(), idx(), ridx() {};

	///
	/// \brief Constructor.
	///
	/// Create the edge (i->j) in the graph.
	/// \param i 	Head of the edge
	/// \param j 	Tail of the edge
	/// \param eij 	Index of the edge
	/// \param eji 	Index of the reversed edge
	edge_id(index i, index j, index eij, index eji) 
		: first(i), second(j), idx(eij), ridx(eji) {};

	///
	/// \brief Indexing operator.
	///
	operator edge_t() const { 
		return edge_t(first,second); 
	}

	static const edge_id NO_EDGE; ///< A static constant for representing missing edges.

	///
	/// \brief Comparison (equal) operator.
	///
	bool operator==(const edge_id& e) const { 
		return e.first==first && e.second==second; 
	};
	
	///
	/// \brief Comparison (not equal) operator.
	///
	bool operator!=(const edge_id& e) const { 
		return !(*this==e); 
	};
	
	///
	/// \brief Comparison (less than) operator.
	///
	bool operator< (const edge_id& e) const { 
		return (first < e.first || (first==e.first && second<e.second)); 
	};
};

typedef std::vector<edge_t> rooted_tree_t;

///
/// \brief The graph structure used by Merlin.
///
/// The graph contains N *nodes* (indexed 0..N-1) and E *edges*, with bidirectional
/// indexing (0..2E-1).
///
class graph {
public:
	typedef edge_id::index index;		///< Basic indexing type for the graph.
	std::vector<set<edge_id> > m_adj;	///< Look up EdgeID info by adj[i][jj] (jj = position of j).
	std::vector<edge_id> m_edges;		///< Look up EdgeID info by edges[eij] (edge index).

protected:
	std::stack<index> m_vvacant;		///< List of available vertex ids.
	std::stack<index> m_evacant;		///< List of available edge ids.

public:
	///
	/// \brief Constructs a graph with *n* vertices and 0 edges.
	/// \param n 	The number of nodes
	///
	explicit graph(size_t n=0) : m_adj(), m_edges(), m_vvacant(), m_evacant() {
		m_adj.resize(n);
	};

	///
	/// \brief Destroys the graph.
	///
	~graph() {
	};

	///
	/// \brief Add a node to the graph.
	/// \return the index of the newly added node.
	///
	index add_node() {			// verify that the adjacency table is large enough
		index use;				// and get an index off the stack if available
		if (m_vvacant.empty()) {
			use = num_nodes();
			m_adj.resize(num_nodes()+1);
		} else {
			use = m_vvacant.top(); m_vvacant.pop();
		}
		return use;
	};

	///
	/// \brief Remove a node from the graph.
	/// \param i 	Index of the node to be removed
	///
	void remove_node(index i) {		// pop_back, or add to vVacant
		m_vvacant.push(i); 	// can we keep track of "real" nNodes and iterate through them?
	};

	///
	/// \brief Return the number of nodes.
	///
	size_t num_nodes() {
		return m_adj.size();
	};

	///
	/// \brief Return the number of edges.
	///
	size_t num_edges() {
		return m_edges.size()/2;
	};

	///
	/// \brief Add an edge to the graph. 

	/// Create the pairs (i,j) and (j,i) in the adjacency list.
	/// \param i 	The head of the edge
	/// \param j 	The tail of the edge
	/// \return the id (index) of the newly added edge.
	///
	const edge_id& add_edge(index i, index j) { // add edges (i,j) and (j,i) to adj
		//std::cout<<"Add edge "<<i<<","<<j<<"\n";
		if (edge(i,j) != edge_id::NO_EDGE)
			return edge(i,j);	// if exists already do nothing

		size_t eij, eji, emax=2*num_edges();	// otherwise get two edge indices
		if (m_evacant.empty()) {
			eij=emax++;
			m_edges.resize(emax);
		} else {
			eij=m_evacant.top();
			m_evacant.pop();
		}
		if (m_evacant.empty()) {
			eji=emax++;
			m_edges.resize(emax);
		} else {
			eji=m_evacant.top();
			m_evacant.pop();
		}

		m_edges[eij] = edge_id(i,j,eij,eji);
		m_adj[i] |= m_edges[eij];
		m_edges[eji] = edge_id(j,i,eji,eij);
		m_adj[j] |= m_edges[eji];

		return m_edges[eij];
	};

	///
	/// \brief Remove an edge from the graph.
	/// \param i 	The head of the edge to be removed
	/// \param j 	The tail of the edge to be removed
	///
	void remove_edge(index i, index j) {	// remove edges (i,j) and (j,i) from adj; add to eVacant
		if (edge(i,j) != edge_id::NO_EDGE) {
			edge_id e = edge(i,j);
			index eij=e.idx, eji = e.ridx;
			m_evacant.push(eij);
			m_evacant.push(eji);// add eij, eji to edge stack
			m_adj[i] /= m_edges[eij];
			m_edges[eij]=edge_id::NO_EDGE;
			m_adj[j] /= m_edges[eji];
			m_edges[eji]=edge_id::NO_EDGE;
		}
	};

	///
	/// \brief Clear the graph by removing all nodes and edges.
	///
	void clear() {
		clear_edges();
		while (!m_vvacant.empty()) m_vvacant.pop();
		m_adj.resize(0);
	};

	///
	/// \brief Remove all edges.
	///
	void clear_edges() {
		m_edges.clear();
		size_t n = num_nodes(); m_adj.clear(); m_adj.resize(n);
		while (!m_evacant.empty()) m_evacant.pop();
	};

	///
	/// \brief Return the edge (if any) between two nodes.
	/// \param i 	The head of the edge
	/// \param j 	The tail of the edge
	/// \return the index of the edge if the edge is present or 
	/// 	the static constant NO_EDGE if the edge is missing.
	///
	const edge_id& edge(index i, index j) const {
		if (i>=m_adj.size())    return edge_id::NO_EDGE;
		set<edge_id>::const_iterator it = m_adj[i].find( edge_id(i,j,0,0) );
		if (it==m_adj[i].end()) return edge_id::NO_EDGE;
		else                   return *it;
	};

	///
	/// \brief Return the neighbors of a node (ie, adjacency list).
	/// \param i 	The index of the node
	/// \return the set of edge ids that are adjacent to the node.
	///
	const set<edge_id>& neighbors(index i) const {
		return m_adj[i];
	};

	///
	/// \brief Return an edge (by its index).
	/// \param eij 	The index of the edge
	/// \return the edge id corresponding to the index.
	const edge_id& edge(index eij) const {
		return m_edges[eij];
	};

	///
	/// \brief Return the edges of the graph.
	/// \return the list of all edge ids in the graph.
	///
	const std::vector<edge_id>& edges() const {
		return m_edges;
	};

	///
	/// \brief Triangulate the graph along an ordering
	/// \param ordering		The variable ordering
	///
	void triangulate(const std::vector<size_t>& ordering);

	///
	/// \brief Create a graph from a set of factor scopes
	/// \param fin		The set of factor scopes
	///
	void init(const vector<variable_set>& fin);

	///
	/// \brief Retrieve the cliques of a triangulated graph
	/// \param ordering		The variable ordering
	///
	std::vector<std::set<size_t> > maximal_cliques(const std::vector<size_t>& ordering);
};

///
/// \brief Output operator.
///
/// Write out the (formatted) content of an edge.
/// \param out 	The output stream to be written
/// \param e 	The edge to be written out
/// \return a reference to the modified output stream containing the edge
///		contents.
///
inline std::ostream& operator<<( std::ostream& out, const edge_t& e) {
	return out << "(" << (int)e.first << "," << (int)e.second <<")";
}

///
/// \brief Output operator.
///
/// Write out the (formatted) content of an edge.
/// \param out 	The output stream to be written
/// \param e 	The edge id to be written out
/// \return a reference to the modified output stream containing the edge
///		contents.
///
inline std::ostream& operator<<( std::ostream& out, const edge_id& e) {
	return out<<(edge_t)e;
}

} // namespace

#endif /* IBM_MERLIN_GRAPH_H_ */
