/*
 * cte.cpp
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


#include "cte.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

using namespace boost;
typedef adjacency_list<vecS, vecS, undirectedS, no_property,
		property<edge_weight_t, int> > BoostGraph;
typedef graph_traits<BoostGraph>::edge_descriptor BoostEdge;
typedef graph_traits<BoostGraph>::vertex_descriptor BoostVertex;
typedef std::pair<int, int> E;

namespace merlin {

// Initialize the clique tree
void cte::init() {

	// Prologue
	if (m_verbose > 0) {
		std::cout << "[CTE] + inference task   : " << m_task << std::endl;
	}

	if (m_query.empty() == false) {
		if (m_verbose > 0) {
			std::cout << "[CTE] + query vars       : ";
			std::copy(m_query.begin(), m_query.end(), std::ostream_iterator<vindex>(std::cout, " "));
			std::cout << std::endl;
		}

		// Add a dummy factor over the query variables
		if (m_query.size() <= MERLIN_MAXSIZE_JOINT_MARGINAL) {
			variable_set scope;
			for (size_t i = 0; i < m_query.size(); ++i) {
				scope |= m_gmo.var(m_query[i]);
			}

			factor q(scope, 1.0);
			m_gmo.add_factor(q);
		}
	}

	if (m_order.size() == 0) { // if we need to construct an elimination ordering
		m_order = m_gmo.order(m_order_method);
	}


	if (m_verbose > 0) {
		std::cout << "[CTE] + ordering method  : " << m_order_method << std::endl;
		std::cout << "[CTE] + elimination      : ";
		std::copy(m_order.begin(), m_order.end(),
				std::ostream_iterator<size_t>(std::cout, " "));
		std::cout << std::endl;

		// Calculate the induced width of the elimination ordering
		size_t wstar = m_gmo.induced_width(m_order);
		std::cout << "[CTE] + induced width    : " << wstar << std::endl;
		std::cout << "[CTE] + exact inference  : Yes" << std::endl;
		std::cout << "[CTE] + ordering time    : " << (timeSystem() - m_start_time) << " seconds" << std::endl;
		std::cout << "[CTE] Building clique tree ... " << std::endl;
	}

	// Get the factors scopes
	size_t n = m_gmo.num_nodes();
	vector<variable_set> fin;
	for (vector<factor>::const_iterator i = m_gmo.get_factors().begin();
			i != m_gmo.get_factors().end(); ++i) {
		fin.push_back((*i).vars());
	}

	// Create the graph
	graph g(n);
	g.init(fin);

	// Triangulate the graph
	g.triangulate(m_order);

	// Find the maximal cliques
	std::vector<std::set<size_t> > clusters = g.maximal_cliques(m_order);

	// Build the clique tree
	build_clique_tree(clusters);

	if (m_verbose > 0) {
		std::cout << "[CTE] Created clique tree with " << m_clusters.size()
				<< " clique factors" << std::endl;
	}

	size_t max_clique_size = 0, max_sepset_size = 0;
	for (size_t i = 0; i < m_edges.size(); ++i) {
		detail::edge* e = m_edges[i];
		max_clique_size = std::max(max_clique_size, e->first->clique.size());
		max_clique_size = std::max(max_clique_size, e->second->clique.size());
		max_sepset_size = std::max(max_sepset_size, e->sepset.size());
	}

	// Initialize beliefs (marginals)
	m_logz = 0;
	m_beliefs.clear();
	m_beliefs.resize(m_gmo.nvar(), factor(1.0));

	// Output the bucket tree statistics
	if (m_verbose > 0) {
		double elapsed = (timeSystem() - m_start_time);
		std::cout << "[CTE] Number of cliques  : " << m_clusters.size() << std::endl;
		std::cout << "[CTE] Number of edges    : " << m_edges.size() << std::endl;
		std::cout << "[CTE] Max clique size    : " << max_clique_size << std::endl;
		std::cout << "[CTE] Max separator size : " << max_sepset_size << std::endl;
		std::cout << "[CTE] Finished initialization in " << elapsed << " seconds" << std::endl;
	}
}

// Build the join tree and initialize the cluster factors
void cte::build_clique_tree(std::vector<std::set<size_t> >& cliques) {

	size_t n = cliques.size();
	std::vector<std::vector<int> > adjacencies;
	adjacencies.resize(n);
	for (std::vector<std::vector<int> >::iterator it = adjacencies.begin();
			it != adjacencies.end(); ++it) {
		it->resize(n);
	}

	// create nodes
	m_clusters.resize(n);
	for (size_t i = 0; i < n; ++i) {
		m_clusters[i].id = i;
		for (std::set<size_t>::iterator j = cliques[i].begin();
				j != cliques[i].end(); ++j) {
			m_clusters[i].clique |= var(*j);
		}
	}

	// fill in the weights
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			if (i != j) {
				detail::node& ni = m_clusters[i];
				detail::node& nj = m_clusters[j];

				std::set<int> sep;
				std::set_intersection(ni.clique.begin(), ni.clique.end(),
						nj.clique.begin(), nj.clique.end(),
						std::inserter(sep, sep.begin()));

				if (sep.empty()) {
					adjacencies[i][j] = (100000);
				} else {
					adjacencies[i][j] = -(sep.size());
				}
			}
		}
	}

	// create a boost weighted graph and run Kruskal on it
	const int num_nodes = n;
	std::vector<E> temp1;
	std::vector<int> temp2;
	for (int i = 0; i < num_nodes - 1; ++i) {
		for (int j = i + 1; j < num_nodes; ++j) {
			if (i != j) {
				temp1.push_back( E(i, j) );
				temp2.push_back( adjacencies[i][j] );
			}
		}
	}

	E* edge_array = new E[temp1.size()];
	int* weights = new int[temp1.size()];
	for (size_t i = 0; i < temp1.size(); ++i) {
		edge_array[i] = temp1[i];
		weights[i] = temp2[i];
	}

	std::size_t num_edges = temp1.size();

	BoostGraph g(edge_array, edge_array + num_edges, weights, num_nodes);
	std::vector<BoostEdge> spanning_tree;

	kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

	if (m_debug) {
		std::cout << "Undirected junction tree (MST):" << std::endl;
		property_map<BoostGraph, edge_weight_t>::type weight = get(edge_weight, g);
		for (std::vector<BoostEdge>::iterator ei = spanning_tree.begin();
				ei != spanning_tree.end(); ++ei) {
			std::cout << source(*ei, g) << " <--> " << target(*ei, g)
					<< " with weight of " << weight[*ei] << "\n";
		}
	}

	// select the root and redirect edges outwards
	m_root = &(m_clusters.back());
	size_t root = m_root->id;
	std::set<size_t> visited;
	std::stack<size_t> dfs;
	dfs.push(root);
	while ( !dfs.empty() ) {

		size_t c = dfs.top(); // current node
		dfs.pop();

		// get all unvisited children and direct them towards n
		for (std::vector<BoostEdge>::iterator ei = spanning_tree.begin();
				ei != spanning_tree.end(); ++ei) {
			size_t src = (size_t) source(*ei, g);
			size_t trg = (size_t) target(*ei, g);
			detail::node* from = NULL, *to = NULL;

			if (src == c &&
				std::find(visited.begin(), visited.end(), trg) == visited.end()) {

				from = &(m_clusters[trg]); // from
				to = &(m_clusters[src]); // to
				dfs.push(trg);
			}

			if (trg == c &&
				std::find(visited.begin(), visited.end(), src) == visited.end()) {

				from = &(m_clusters[src]); // from
				to = &(m_clusters[trg]); // to
				dfs.push(src);
			}

			// create the directed edge
			if (from != NULL && to != NULL) {
				detail::edge* e = new detail::edge(from, to);
				m_edges.push_back(e);
				from->edges.push_back(e);
				to->edges.push_back(e);
				to->children.push_back(from);
				from->parent = to;
			}
		}

		visited.insert(c);
	}

	// find the message order (schedule) by a bfs of the junction tree
	std::queue<detail::node*> bfs;
	bfs.push(m_root);
	while (!bfs.empty()) {
		detail::node* c = bfs.front();
		bfs.pop();
		vector<detail::edge*>& elist = c->edges;
		for (vector<detail::edge*>::iterator it = elist.begin();
				it != elist.end(); ++it) {
			detail::edge* e = (*it);
			if (e->second->id == c->id) {
				m_messages.push_back(e);
			}
		}

		for (vector<detail::node*>::iterator it = c->children.begin();
				it != c->children.end(); ++it) {
			bfs.push(*it);
		}
	}

	// safety checks
	assert(m_messages.size() == m_edges.size());

	// reverse the order to start from the leaves
	std::reverse(m_messages.begin(), m_messages.end());

	// clean up
	delete[] weights;
	delete[] edge_array;

	if (m_debug) {
		std::ofstream fout("kruskal.dot");
		fout << "digraph JT {\n" << " rankdir=LR\n" << " size=\"3,3\"\n"
				<< " ratio=\"filled\"\n" << " edge[style=\"bold\"]\n"
				<< " node[shape=\"circle\"]\n";
		// nodes
		for (size_t i = 0; i < m_clusters.size(); ++i) {
			std::stringstream label;
			label << m_clusters[i].id << ": " << m_clusters[i].clique;
			fout << "node" << m_clusters[i].id
				<< "[ label = \"" << label.str() << "\"];\n";
		}

		// edges
		for (size_t i = 0; i < m_edges.size(); ++i) {
			detail::edge* e = m_edges[i];
			fout << "node" << e->first->id << " -> " << "node" << e->second->id << ";\n";
		}

		fout << "}\n";
		fout.close();
	}

	// Map variables to clusters
	m_var2clique.resize(m_gmo.nvar(), -1);
	for (size_t i = 0; i < m_clusters.size(); ++i) {
		detail::node& cl = m_clusters[i];
		variable_set& vars = cl.clique;
		for (variable_set::const_iterator j = vars.begin();
				j != vars.end(); ++j) {
			size_t v = (*j).label();
			if (m_var2clique[v] < 0) {
				m_var2clique[v] = cl.id;
			}
		}
	}

	// Allocate original functions to clusters and compute clique factors
	size_t idx = 0;
	for (vector<factor>::const_iterator fi = m_gmo.get_factors().begin();
			fi != m_gmo.get_factors().end(); ++fi, ++idx) {

		const factor& f = (*fi);
		for (size_t ci = 0; ci != m_clusters.size(); ++ci) {
			detail::node& cl = m_clusters[ci];
			if ((f.vars() & cl.clique) == f.vars()) {
				cl.originals.push_back(idx);
				cl.theta *= f;
				break;
			}
		}
	}

	if (m_debug) {
		std::cout << "Initial cluster factors:" << std::endl;
		for (size_t ci = 0; ci != m_clusters.size(); ++ci) {
			detail::node& cl = m_clusters[ci];
			std::cout << " " << ci << " -> " << cl.theta << std::endl;
		}
	}
}

// Update the collection of factors, recompute the clique factors and reset the
// messages. We assume there is a one-to-one mapping between the old factors
// and the new factors. This method is to be used only for EM param learning.
void cte::reinit(const std::vector<factor>& factors) {

	// Set the new factors
	for (size_t i = 0; i < factors.size(); ++i) {
		m_gmo.set_factor(i, factors[i]);
	}

	// Recompute the clique factors
	for (size_t ci = 0; ci != m_clusters.size(); ++ci) {
		detail::node& cl = m_clusters[ci];
		cl.theta = factor(1.0);
		for (size_t j = 0; j < cl.originals.size(); ++j) {
			const factor& f = m_gmo.get_factor(cl.originals[j]);
			cl.theta *= f;
		}
	}

	// Reset the edge messages
	for (vector<detail::edge*>::iterator it = m_edges.begin();
			it != m_edges.end(); ++it) {
		(*it)->reset();
	}

	// Reset the log partition function and beliefs
	m_logz = 0;
	m_beliefs.clear();
	m_beliefs.resize(m_gmo.nvar(), factor(1.0));

	if (m_debug) {
		std::cout << "Initial cluster factors:" << std::endl;
		for (size_t ci = 0; ci != m_clusters.size(); ++ci) {
			detail::node& cl = m_clusters[ci];
			std::cout << " " << ci << " -> " << cl.theta << std::endl;
		}
	}
}

// Forward message propagation (from leaves to root)
void cte::forward() {
	double timestamp = timeSystem();
	for (vector<detail::edge*>::iterator i = m_messages.begin();
			i != m_messages.end(); ++i) {
		detail::edge* e = (*i);
		e->messageFwd();
	}

	if (m_verbose > 0) {
		std::cout << "[CTE] Finished forward pass in " << (timeSystem() - timestamp)
			<< " seconds" << std::endl;
	}
}

// Forward message propagation with evidence (from leaves to root)
void cte::forward(std::vector<int>& evidence) {
	for (vector<detail::edge*>::iterator i = m_messages.begin();
			i != m_messages.end(); ++i) {
		detail::edge* e = (*i);
		e->messageFwd(evidence);

		if (m_debug) {
			std::cout << " -> forward msg from "
				<< e->first->id << " to " << e->second->id
				<< ": " << e->fwd << std::endl;
		}
	}
}

// Forward message on an edge
void detail::edge::messageFwd() {

	// collect the original functions
	factor F = first->theta;

	// collect the incoming messages (should already be projected)
	for (size_t i = 0; i < first->edges.size(); ++i) {
		detail::edge* e = first->edges[i];
		if (e == this) {
			continue; // skip current edge
		}
		if (e->first->id == first->id) {
			F *= e->bwd; //e->getMessage2();
		} else if (e->second->id == first->id) {

			F *= e->fwd; //e->getMessage1();
		}
	}

	// marginalize the eliminator
	variable_set elim = (first->clique - sepset);
	if (elim.size() > 0) {
		F = F.sum(elim);
	}

	// update the forward message
	this->fwd = F;
}

// Forward message computed with evidence on an edge
void detail::edge::messageFwd(std::vector<int>& evidence) {

	// collect the original functions
	factor F = first->theta;
	F = F.condition(evidence);

	// collect the incoming messages (should already be projected)
	for (size_t i = 0; i < first->edges.size(); ++i) {
		detail::edge* e = first->edges[i];
		if (e == this) {
			continue; // skip current edge
		}
		if (e->first->id == first->id) {
			F *= e->bwd.condition(evidence);
		} else if (e->second->id == first->id) {

			F *= e->fwd.condition(evidence);
		}
	}

	// Find the evidence variables in the scope of the cluster
	variable_set evid;
	for (variable_set::const_iterator i = first->clique.begin();
			i != first->clique.end(); ++i) {
		int v = i->label();
		if (evidence[v] >= 0) {
			evid |= *i;
		}
	}

	// Marginalize the non-evidence eliminator
	//F = F.condition(evidence);
	variable_set elim = (first->clique - sepset);
	elim -= evid;
	if (elim.size() > 0) {
		F = F.sum(elim);
	}

	// update the forward message
	this->fwd = F;

}

// Backward message propagation (from root to leaves)
void cte::backward() {
	double timestamp = timeSystem();
	for (vector<detail::edge*>::reverse_iterator ri = m_messages.rbegin();
			ri != m_messages.rend(); ++ri) {
		detail::edge* e = (*ri);
		e->messageBwd();
	}

	if (m_verbose > 0) {
		std::cout << "[CTE] Finished backward pass in " << (timeSystem() - timestamp)
			<< " seconds" << std::endl;
	}
}

// Backward message propagation with evidence (from root to leaves)
void cte::backward(std::vector<int>& evidence) {
	for (vector<detail::edge*>::reverse_iterator ri = m_messages.rbegin();
			ri != m_messages.rend(); ++ri) {
		detail::edge* e = (*ri);
		e->messageBwd(evidence);

		if (m_debug) {
			std::cout << " <- backward msg from "
				<< e->second->id << " to " << e->first->id
				<< ": " << e->bwd << std::endl;
		}
	}
}

// Backward message on an edge
void detail::edge::messageBwd() {

	// collect the original functions
	factor F = second->theta;

	// collect the incoming messages
	for (size_t i = 0; i < second->edges.size(); ++i) {
		detail::edge* e = second->edges[i];
		if (e == this) {
			continue; // skip current edge
		}
		if (e->first->id == second->id) {
			F *= e->bwd; //e->getMessage2();
		} else if (e->second->id == second->id) {
			F *= e->fwd; //e->getMessage1();
		}
	}

	// marginalize the eliminator
	variable_set elim = (second->clique - sepset);
	if (elim.size() > 0) {
		F = F.sum(elim);
	}

	this->bwd = F;
}

// Backward message computed with evidence on an edge
void detail::edge::messageBwd(std::vector<int>& evidence) {

	// collect the original functions
	factor F = second->theta;
	F = F.condition(evidence);

	// collect the incoming messages
	for (size_t i = 0; i < second->edges.size(); ++i) {
		detail::edge* e = second->edges[i];
		if (e == this) {
			continue; // skip current edge
		}
		if (e->first->id == second->id) {
			F *= e->bwd.condition(evidence); //e->getMessage2();
		} else if (e->second->id == second->id) {
			F *= e->fwd.condition(evidence); //e->getMessage1();
		}
	}

	// find the evidence variables in the scope of the cluster
	variable_set evid;
	for (variable_set::const_iterator i = second->clique.begin();
			i != second->clique.end(); ++i) {
		int v = i->label();
		if (evidence[v] >= 0) {
			evid |= *i;
		}
	}

	// marginalize the eliminator
	variable_set elim = (second->clique - sepset);
	elim -= evid;
	if (elim.size() > 0) {
		F = F.sum(elim);
	}

	this->bwd = F;

}

// Clique tree calibration
void cte::calibrate() {

	forward();
	backward();
}

// Update the beliefs (marginals or max-marginals)
void cte::update() {

	// Compute the clique beliefs
	for (size_t i = 0; i < m_clusters.size(); ++i) {
		detail::node& cl = m_clusters[i];
		cl.belief = cl.theta;
		for (vector<detail::edge*>::iterator ei = cl.edges.begin();
				ei != cl.edges.end(); ++ei) {
			detail::edge* e = (*ei);
			if (e->second->id == cl.id) {
				cl.belief *= e->fwd;
			}
			if (e->first->id == cl.id) {
				cl.belief *= e->bwd;
			}
		}
	}

	// Compute the log partition function
	m_logz = std::log(m_root->belief.sum());

	// Compute the marginals
	for (size_t v = 0; v < m_var2clique.size(); ++v) {
		size_t ci = m_var2clique[v];
		detail::node& cl = m_clusters[ci];
		variable VX = m_gmo.var(v);

		m_beliefs[v] = marg(cl.belief, VX);
		m_beliefs[v].normalize(); // normalize
	}
}

// Compute the joint marginal of a set of variables
void cte::joint_marginal(const variable_set& scope) {

	// Find the shallowest cliques that contain the scope
	std::vector<int> nodes;
	variable_set temp = scope;
	while (temp.size() > 0) {
		int max_id = -1, max_score = -1;
		std::queue<detail::node*> bfs;
		bfs.push(m_root);
		while (!bfs.empty()) {
			detail::node* c = bfs.front();
			bfs.pop();

			int score = (c->clique & temp).size();
			if (score > 0 && score > max_score) {
				max_score = score;
				max_id = c->id;
			}

			for (vector<detail::node*>::iterator it = c->children.begin();
					it != c->children.end(); ++it) {
				bfs.push(*it);
			}
		}

		assert(max_id >= 0);
		nodes.push_back(max_id);
		temp = temp - (temp & m_clusters[max_id].clique);
	}

	if (m_debug) {
		std::cout << "[DEBUG] Found the shallowest nodes: ";
		std::copy(nodes.begin(), nodes.end(),
				std::ostream_iterator<int>(std::cout, " "));
		std::cout << std::endl;
	}

	// Check if the scope is included in one clique
	if (nodes.size() == 1) { // one clique
		m_marginal = marg(m_clusters[nodes.at(0)].belief, scope);
		m_marginal.normalize();
	} else { // multiple cliques

		// Collect all relevant factors
		vector<factor> factors;
		factors.push_back(m_root->belief);
		for (size_t i = 0; i < nodes.size(); ++i) {
			detail::node* c = &(m_clusters.at(nodes[i]));
			while (c != NULL) {
				for (vector<detail::edge*>::iterator ei = c->edges.begin();
						ei != c->edges.end(); ++ei) {
					detail::edge* e = (*ei);
					if (e->first->id == c->id) {
						factor f = c->belief;
						f /= e->fwd;
						factors.push_back(f);
						break;
					}
				}

				c = c->parent;
				if (c == m_root) {
					break;
				}
			}
		}

		variable_set all_vars, elim_vars;
		for (vector<factor>::iterator fi = factors.begin();
				fi != factors.end(); ++fi) {
			all_vars |= fi->vars();
		}
		elim_vars = (all_vars - scope);

		// Run variables elimination over these factors
		variable_order_t order;
		for (size_t i = 0; i < m_order.size(); ++i) {
			size_t v = m_order[i];
			if (elim_vars.contains(m_gmo.var(v))) {
				order.push_back(v);
			}
		}

		if (m_debug) {
			std::cout << "[DEBUG] All vars: " << all_vars << std::endl;
			std::cout << "[DEBUG] Scope: " << scope << std::endl;
			std::cout << "[DEBUG] Elim: " << elim_vars << std::endl;
			std::cout << "[DEBUG] Factors: " << factors.size() << std::endl;
			std::cout << "[DEBUG] Elim order: ";
			std::copy(order.begin(), order.end(),
					std::ostream_iterator<size_t>(std::cout, " "));
			std::cout << std::endl;
		}

		// Run variable elimination
		for (size_t i = 0; i < order.size(); ++i) {
			size_t v = order[i];
			variable VX = m_gmo.var(v);

			// Collect all factors mentioning VX
			factor f(1.0);
			for (vector<factor>::iterator it = factors.begin();
					it != factors.end();) {
				if (it->variables().contains(VX)) {
					f *= (*it);
					it = factors.erase(it);
				} else {
					++it;
				}
			}

			// Eliminate the variable
			f = elim(f, VX);
			factors.push_back(f);
		}

		// Compute the joint marginal
		m_marginal = factor(1.0);
		for (vector<factor>::iterator it = factors.begin();
				it != factors.end(); ++it) {
			m_marginal *= (*it);
			m_marginal.normalize();
		}
	}
}

// The scope must be included in one cluster (Bayes nets only)
void cte::joint_marginal(const variable_set& scope, std::vector<int>& evidence) {

	// Get the cluster that contains the scope
	bool found = false;
	size_t j = 0;
	for (size_t i = 0; i < m_clusters.size(); ++i) {
		detail::node& cl = m_clusters[i];
		variable_set temp = (cl.clique & scope);
		if (temp == scope) {
			found = true;
			j = i;
			break;
		}
	}

	// Safety check
	assert(found);

	// Compute the joint marginal (belief) subject to evidence
	m_marginal = factor(scope, 0.0);
	detail::node& n = m_clusters[j];
	n.belief = n.theta.condition(evidence);
	for (vector<detail::edge*>::iterator ei = n.edges.begin();
			ei != n.edges.end(); ++ei) {
		detail::edge* e = (*ei);
		if (e->second->id == n.id) {
			n.belief *= e->fwd;
		}
		if (e->first->id == n.id) {
			n.belief *= e->bwd;
		}
	}

	// Update the belief entries subject to evidence
	variable_set v_rem;
	for (variable_set::const_iterator vi = n.clique.begin();
			vi != n.clique.end(); ++vi) {
		if (!scope.contains(*vi) && evidence[vi->label()] < 0) {
			v_rem |= *vi;
		}
	}
	n.belief = n.belief.sum(v_rem);

	if (m_debug) {
		std::cout << "[DEBUG] Joint marginal scope: " << m_marginal.vars() << std::endl;
		std::cout << "[DEBUG] Actual belief scope:  " << n.belief << std::endl;
	}

	// Fill in all configurations of the marginal (including the evidence)
	index_config cv1(m_marginal.vars(), true);
	config_index cv2(n.belief.vars(), true);
	for (size_t i = 0; i < m_marginal.num_states(); ++i) {
		std::map<size_t, size_t> config = cv1.convert(i);
		if (is_compatible(config, evidence)) {
			size_t j = cv2.convert(config);
			double v = n.belief.get(j);
			m_marginal.set(i, v);
		} else {
			m_marginal.set(i, 0.0);
		}
	}

	// Normalize by P(evidence)
	m_marginal /= std::exp(m_logz);

	if (m_debug) {
		std::cout << "[DEBUG] Joint marginal: " << m_marginal << std::endl;
	}
}

// Check if a variable configuration is compatible with the evidence vector
bool cte::is_compatible(std::map<size_t, size_t>& config,
		std::vector<int>& evidence) {
	bool ok = true;
	for (std::map<size_t, size_t>::iterator mi = config.begin();
			mi != config.end(); ++mi) {
		size_t var = mi->first;
		size_t val = mi->second;
		if (evidence[var] < 0) {
			continue; // skip non-evidence variable
		} else {
			if (val != evidence[var]) {
				ok = false;
				break;
			}
		}
	}

	return ok;
}

// Run the clique-tree algorithm with evidence (for EM learning)
// Returns *true* if P(evidence) > 0, and *false* otherwise
bool cte::propagate_evidence(std::vector<int>& evidence) {

	// Save the current evidence
	if (m_debug) {
		std::cout << "[CTE] Propagate evidence: ";
		std::copy(evidence.begin(), evidence.end(), std::ostream_iterator<int>(std::cout, " "));
		std::cout << std::endl;
	}

	m_evidence = evidence;
	forward(evidence);
	backward(evidence);

	// Compute the posterior marginals for all variables (evidence, non-evidence)
	for (size_t i = 0; i < m_clusters.size(); ++i) {
		detail::node& cl = m_clusters[i];
		cl.belief = cl.theta.condition(evidence);
		for (vector<detail::edge*>::iterator ei = cl.edges.begin();
				ei != cl.edges.end(); ++ei) {
			detail::edge* e = (*ei);
			if (e->second->id == cl.id) {
				cl.belief *= e->fwd;
			}
			if (e->first->id == cl.id) {
				cl.belief *= e->bwd;
			}
		}
	}

	// Compute the root clique belief (to get log partition function)
	m_root->belief = m_root->theta.condition(evidence);
	for (vector<detail::edge*>::iterator ei = m_root->edges.begin();
		ei != m_root->edges.end(); ++ei) {
		detail::edge* e = (*ei);
		if (e->second->id == m_root->id) {
			m_root->belief *= e->fwd;
		}
		if (e->first->id == m_root->id) {
			m_root->belief *= e->bwd;
		}
	}

	// Compute the probability of evidence
	bool result = true;
	factor::value pe = m_root->belief.sum();
	if (pe == 0.0) {
		result = false;
	}

	// Update the log partition function
	m_logz = std::log(m_root->belief.sum());

	// Compute the marginals
	for (size_t v = 0; v < m_var2clique.size(); ++v) {
		size_t ci = m_var2clique[v];
		detail::node& cl = m_clusters[ci];
		variable VX = m_gmo.var(v);

		if (evidence[v] >= 0) {
			factor tmp(variable_set(VX), 0.0);
			tmp.set(evidence[v], 1.0);
			m_beliefs[v] = tmp;
		} else {
			m_beliefs[v] = marg(cl.belief, VX);
			m_beliefs[v].normalize(); // normalize
		}
	}

	if (m_debug) {
		std::cout << "[CTE] Finished propagating evidence with logZ = "
				<< m_logz << " (" << std::exp(m_logz) << ")" << std::endl;
		std::cout << "[CTE] Posterior marginals:" << std::endl;
		for (size_t i = 0; i < m_beliefs.size(); ++i) {
			std::cout << " " << m_beliefs[i] << std::endl;
		}
	}

	return result;
}

/// Run the clique-tree elimination algorithm.
void cte::run() {

	// Start the timer and store it
	m_start_time = timeSystem();

	init();
	calibrate();
	update();
	if (!m_query.empty()) {
		variable_set scope;
		for (size_t i = 0; i < m_query.size(); ++i) {
			scope |= m_gmo.var(m_query[i]);
		}

		joint_marginal(scope);
	}

	// Output solution (UAI output format)
	if (m_verbose > 0) {
		std::cout << "[CTE] Finished in " << (timeSystem() - m_start_time)
			<< " seconds" << std::endl;
	}

	// Compute the joint marginal given in query

	switch (m_task) {
	case Task::PR:	// partition function
		{
			std::cout << "PR" << std::endl;
			std::cout << std::fixed << std::setprecision(MERLIN_PRECISION)
				<< m_logz << " (" << std::scientific
				<< std::setprecision(MERLIN_PRECISION)
				<< std::exp(m_logz) << ")" << std::endl;

			if (isinf(m_logz)) {
				std::cout << "STATUS" << std::endl;
				std::cout << "false: Inconsistent evidence or underflow" << std::endl;
			} else {
				std::cout << "STATUS" << std::endl;
				std::cout << "true: Consistent evidence" << std::endl;
			}

			break;
		}
	case Task::MAR:	// marginals
		{
			std::cout << "PR" << std::endl;
			std::cout << std::fixed << std::setprecision(MERLIN_PRECISION)
				<< m_logz << " (" << std::scientific
				<< std::setprecision(MERLIN_PRECISION)
				<< std::exp(m_logz) << ")" << std::endl;

			if (isinf(m_logz)) {
				std::cout << "STATUS" << std::endl;
				std::cout << "false: Inconsistent evidence or underflow" << std::endl;
			} else {
				std::cout << "STATUS" << std::endl;
				std::cout << "true: Consistent evidence" << std::endl;
			}

			std::cout << "MAR" << std::endl;
			std::cout << m_gmo.nvar();
			for (vindex v = 0; v < m_gmo.nvar(); ++v) {
				variable VX = m_gmo.var(v);
				std::cout << " " << VX.states();
				for (size_t k = 0; k < VX.states(); ++k) {
					std::cout << " " << std::fixed
						<< std::setprecision(MERLIN_PRECISION)
						<< belief(VX)[k];
				}
			}
			std::cout << std::endl;
			if (!m_query.empty()) {
				variable_set scope = m_marginal.variables();
				std::cout << "JOINT_MAR : " << scope << std::endl;
				std::vector<size_t> dims(scope.size());
				size_t j = 0;
				for (variable_set::const_iterator si = scope.begin(); si != scope.end(); ++si) {
					dims[j++] = (*si).states();
				}

				for (size_t index = 0; index < m_marginal.numel(); ++index) {
					std::vector<size_t> I(dims.size(), 0); // configuration of source variables
					size_t i = index;
					for (size_t v = 0; v < dims.size(); ++v) {
						I[v] = i % dims[v];
						i -= I[v];
						i /= dims[v];
					}

					for (size_t j = 0; j < I.size(); ++j) {
						std::cout << I[j] << " ";
					}
					std::cout << ": " << std::fixed
							<< std::setprecision(MERLIN_PRECISION)
							<< m_marginal[index]
							<< std::endl;
				}
			}

			break;
		}
	default:
		break;
	}
}

/// Write the solution to the output stream
void cte::write_solution(std::ostream& out, const std::map<size_t, size_t>& evidence,
		const std::map<size_t, size_t>& old2new, const graphical_model& orig,
		const std::set<size_t>& dummies, int output_format) {

	if (output_format == MERLIN_OUTPUT_JSON) {
		out << "{";
		out << " \"algorithm\" : \"cte\", ";
		switch (m_task) {
		case Task::PR:
			{
				double val = m_logz + std::log(orig.get_global_const());
				double prob = std::exp(val);

				out << " \"task\" : \"PR\", ";
				out << " \"value\" : " << std::fixed
					<< std::setprecision(MERLIN_PRECISION)
					<< (m_logz + std::log(orig.get_global_const())) << ", ";
				if (prob == 0.0) { // probability of evidence is 0
					out << " \"status\" : \"false\", ";
					out << " \"message\" : \"Inconsistent evidence or underflow\" ";
				} else {
					out << " \"status\" : \"true\", ";
					out << " \"message\" : \"Consistent evidence\" ";
				}

				break;
			}
		case Task::MAR:
			{
				double val = m_logz + std::log(orig.get_global_const());
				double prob = std::exp(val);

				out << " \"task\" : \"MAR\", ";
				out << " \"value\" : " << std::fixed
					<< std::setprecision(MERLIN_PRECISION)
					<< (m_logz + std::log(orig.get_global_const())) << ", ";

				if (prob == 0.0) { // probability of evidence is 0
					out << " \"status\" : \"false\", ";
					out << " \"message\" : \"Inconsistent evidence or underflow\", ";
					out << " \"marginals\" : [] ";
				} else {
					out << " \"status\" : \"true\", ";
					out << " \"message\" : \"Consistent evidence\", ";
					out << " \"marginals\" : [ ";

					for (vindex i = 0; i < orig.nvar(); ++i) {
						if (dummies.find(i) != dummies.end()) {
							continue; // skip dummy variables from output
						}

						variable v = orig.var(i);
						out << "{";
						out << " \"variable\" : " << v.label() << ", ";
						out << " \"states\" : " << v.states() << ", ";
						out << " \"probabilities\" : [";
						try { // evidence variable
							size_t val = evidence.at(i);
							for (size_t k = 0; k < v.states(); ++k) {
								out << std::fixed
									<< std::setprecision(MERLIN_PRECISION)
									<< (k == val ? 1.0 : 0.0);
								if (k != v.states() - 1) {
									out << ", ";
								}
							}
							out << "] ";
						} catch(std::out_of_range& e) { // non-evidence variable
							vindex vx = old2new.at(i);
							variable VX = var(vx);
							for (size_t k = 0; k < VX.states(); ++k) {
								out << std::fixed
									<< std::setprecision(MERLIN_PRECISION)
									<< belief(VX)[k];
								if (k != VX.states() - 1) {
									out << ", ";
								}
							}
							out << "] ";
						}
						out << "}";
						if (i != orig.nvar() - 1) {
							out << ", ";
						}
					} // end for
					out << "] ";

					if (!m_query.empty()) {
						out << ", ";
						out << "\"joint_marginal\" : {";
						variable_set scope = m_marginal.variables();
						out << "\"scope\" : [";
						for (size_t i = 0; i < m_query.size(); ++i) {
							out << m_query[i];
							if (i < m_query.size() - 1) {
								out << ",";
							}
						}
						out << "], " << std::endl;
						out << "\"probabilities\" : [";

						std::vector<size_t> dims(scope.size());
						size_t j = 0;
						for (variable_set::const_iterator si = scope.begin();
								si != scope.end(); ++si) {
							dims[j++] = (*si).states();
						}

						size_t N = m_marginal.numel();
						for (size_t index = 0; index < N; ++index) {
							std::vector<size_t> I(dims.size(), 0); // configuration of source variables
							size_t i = index;
							for (size_t v = 0; v < dims.size(); ++v) {
								I[v] = i % dims[v];
								i -= I[v];
								i /= dims[v];
							}

							out << "{\"config\" : [";
							for (size_t j = 0; j < I.size(); ++j) {
								out << I[j];
								if (j < I.size() - 1) {
									out << ",";
								}
							}
							out << "], ";
							out << "\"value\" : " << std::fixed
									<< std::setprecision(MERLIN_PRECISION)
									<< m_marginal[index]
									<< "}";
							if (index < N - 1) {
								out << ", ";
							}
						}
						out << "]}";
					}
				}

				break;
			}
		default:
			break;
		}

		out << "}";
	} else if (output_format == MERLIN_OUTPUT_UAI) {
		switch (m_task) {
		case Task::PR:
		case Task::MAR:
			{
				double val = m_logz + std::log(orig.get_global_const());
				double prob = std::exp(val);

				out << "PR" << std::endl;
				out << std::fixed << std::setprecision(MERLIN_PRECISION)
					<< (m_logz + std::log(orig.get_global_const())) << " ("
					<< std::scientific << std::setprecision(MERLIN_PRECISION)
					<< std::exp(m_logz + std::log(orig.get_global_const()))
					<< ")" << std::endl;

				if (prob == 0.0) {
					out << "STATUS" << std::endl;
					out << "false: Inconsistent evidence or underflow" << std::endl;
				} else {
					out << "STATUS" << std::endl;
					out << "true: Consistent evidence" << std::endl;
				}

				out << "MAR" << std::endl;
				out << orig.nvar() - dummies.size();
				for (vindex i = 0; i < orig.nvar(); ++i) {
					if (dummies.find(i) != dummies.end()) {
						continue; // skip dummy variables from output
					}

					variable v = orig.var(i);
					try { // evidence variable
						size_t val = evidence.at(i);
						out << " " << v.states();
						for (size_t k = 0; k < v.states(); ++k) {
							out << " " << std::fixed
								<< std::setprecision(MERLIN_PRECISION)
								<< (k == val ? 1.0 : 0.0);
						}
					} catch(std::out_of_range& e) { // non-evidence variable
						vindex vx = old2new.at(i);
						variable VX = var(vx);
						out << " " << VX.states();
						for (size_t k = 0; k < VX.states(); ++k) {
							out << " " << std::fixed
								<< std::setprecision(MERLIN_PRECISION)
								<< belief(VX)[k];
						}
					}
				} // end for
				out << std::endl;

				if (!m_query.empty()) {
					variable_set scope = m_marginal.variables();
					out << "JOINT_MAR : [";
					for (size_t i = 0; i < m_query.size(); ++i) {
						out << m_query[i];
						if (i < m_query.size() - 1) {
							out << ",";
						}
					}
					out << "]" << std::endl;

					std::vector<size_t> dims(scope.size());
					size_t j = 0;
					for (variable_set::const_iterator si = scope.begin();
							si != scope.end(); ++si) {
						dims[j++] = (*si).states();
					}

					for (size_t index = 0; index < m_marginal.numel(); ++index) {
						std::vector<size_t> I(dims.size(), 0); // configuration of source variables
						size_t i = index;
						for (size_t v = 0; v < dims.size(); ++v) {
							I[v] = i % dims[v];
							i -= I[v];
							i /= dims[v];
						}

						for (size_t j = 0; j < I.size(); ++j) {
							out << I[j] << " ";
						}
						out << ": " << std::fixed
								<< std::setprecision(MERLIN_PRECISION)
								<< m_marginal[index]
								<< std::endl;
					}
				}

				break;
			}
		default:
			break;
		}
	} else {
		std::string err_msg("Unknown output format.");
		throw std::runtime_error(err_msg);
	}
}


} // namespace

