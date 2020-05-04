/*
 * graphical_model.h
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

/// \file graphical_model.h
/// \brief A probabilistic graphical model
/// \author Radu Marinescu radu.marinescu@ie.ibm.com


#ifndef IBM_MERLIN_GRAPHICAL_MODEL_H_
#define IBM_MERLIN_GRAPHICAL_MODEL_H_

#include "enum.h"
#include "factor.h"
#include "graph.h"


namespace merlin {

///
/// \brief Graphical model base class.
///
/// Simplest form of a *graphical model*, namely just a collection of 
/// variables and factors (defined on subsets of variables).
/// Internally, factors and variables are mainly referenced by integer indices, with
///  0 <= f < nFactors() and 0 <= v < nvar().  However, many interface functions are called
///  using a variable object, Var(label,dim).  To convert, var(v) gives the vth Var object
///  and the internal function _vindex(V) gives the index corresponding to variable object V.
///
class graphical_model: public graph {
public:

	// Useful typedefs:

	typedef size_t findex;			///< Factor index
	typedef size_t vindex;			///< Variable index
	typedef set<findex> flist; 		///< Collection of factor indices

	// Constructors, copy and assignment:

	///
	/// \brief Creates an empty graphical model.
	///
	graphical_model() :
		m_factors(), m_markov(true), m_vadj(), m_dims() {
	};

	///
	/// \brief Creates a graphical model by copying an existing one.
	/// \param gm 	An object of the same type
	///
	graphical_model(const graphical_model& gm) :
			graph((graph&) gm), m_factors(gm.m_factors), m_markov(gm.m_markov),
			m_vadj(gm.m_vadj), m_dims(gm.m_dims) {
	};

	///
	/// \brief Assignment operator (deep copy).
	/// \param gm 	The graphical model to be copied from
	/// \return a reference to the modified object containing the copied content.
	///
	graphical_model& operator=(const graphical_model& gm) {
		graph::operator=((graph&) gm);			// copy graph elements over
		m_factors = gm.m_factors;				// copy list of factors
		m_vadj = gm.m_vadj;						// variable dependencies
		m_dims = gm.m_dims;						// and variable dimensions over
		m_markov = gm.m_markov;					// and format
		return *this;
	};

	///
	/// \brief Clone the graphical model.
	/// \return the pointer to the new object representing the copy of 
	/// the current model.
	///
	graphical_model* clone() {
		graphical_model* gm = new graphical_model(*this);
		return gm;
	};

	///
	/// \brief Constructor from a list of factors.
	/// \param fs 	The list of factors
	///
	graphical_model(std::vector<factor> fs) :
			graph(), m_factors(fs), m_vadj(), m_dims() {
		fixup();
	};

	///
	/// \brief Creates a graphical model from input iterators.
	/// \param first 	The iterator to beginning
	/// \param last 	The iterator to end
	///
	template<class InputIterator>
	graphical_model(InputIterator first, InputIterator last) :
			m_factors(first, last), m_vadj(), m_dims() {
		fixup();
	}

	///
	/// \brief Destroys the graphical model.
	///
	virtual ~graphical_model() {
	};

	///
	/// \brief Read the graphical model from a file in the UAI format.
	///
	///	For details on the UAI file format see the main documentation of
	/// the library. In this format, the factor tables are assumed to be
	/// represented using the "least significant bit", namely the last variable
	/// in the factor scope changes the fastest. Moreover, the internal
	/// representation of the factor assumes that the scope is ordered
	/// lexicographically.
	///
	/// \param file_name 	The full path to the file
	///
	void read(std::istream& is, bool positive_mode = false) {

		// Read the header
		size_t nvar, ncliques, csize, v, nval;
		std::string st;
		is >> st;
		if ( st.compare("MARKOV") == 0 ) {
			m_markov = true;
		} else if ( st.compare("BAYES") == 0 ) {
			m_markov = false;
		} else {
			std::string err_msg("Merlin only supports the UAI Markov or Bayes file format.");
			throw std::runtime_error(err_msg);
		}

		// Read the number of variables and their domains
		is >> nvar;
		std::vector<size_t> dims(nvar);
		for (size_t i = 0; i < nvar; i++)
			is >> dims[i];

		// Read the number of factors and their scopes (scope is a variable_set)
		is >> ncliques;
		std::vector<std::vector<variable> > cliques(ncliques);
		std::vector<variable_set> sets(ncliques);
		for (size_t i = 0; i < ncliques; i++) {
			is >> csize;
			cliques[i].reserve(csize);
			for (size_t j = 0; j < csize; j++) {
				is >> v;
				variable V(v, dims[v]);
				cliques[i].push_back(V);
				sets[i] |= V;
			}
		}

		// Read the factor tables (ensure conversion to ordered scopes)
		double fval;
		std::vector<factor> factors(ncliques);
		for (size_t i = 0; i < ncliques; i++) {
			is >> nval;
			assert(nval == sets[i].num_states());
			factors[i] = factor(sets[i], 0.0); // preallocate memory
			if (m_markov == false) { // for Bayes nets, last variable is the child
				factors[i].set_child(cliques[i].back().label());
			}
			convert_index ci(cliques[i], false, true); // convert from source order (littleEndian) to target order (bigEndian)
			for (size_t j = 0; j < nval; j++) {
				size_t k = ci.convert(j);	// get the index in the factor table
				is >> fval; // read the factor value;
				if (positive_mode) { // force positive values (> 0) if enabled
					if (fval == 0.0) {
						fval += MERLIN_EPSILON; // adjust slightly for numerical stability
					} else if (fval == 1.0) {
						fval -= MERLIN_EPSILON; // adjust slightly for numerical stability
					}
				}

				factors[i][k] = fval; // save the factor value into the table
			}
		}

		m_factors = factors;
		fixup();
	}

	///
	/// \brief Write the graphical model to a file using the UAI format.
	///
	/// For details on the UAI file format see the main documentation of
	/// the library. In this format, the factor tables are assumed to be
	/// represented using the "least significant bit", namely the last variable
	/// in the factor scope changes the fastest. Moreover, the internal
	/// representation of the factor assumes that the scope is ordered
	/// lexicographically.
	///
	void write(std::ostream& os) {

		// Write the header
		os << "MARKOV" << std::endl;
		os << nvar() << std::endl;
		for (size_t i = 0; i < m_dims.size(); ++i) {
			os << m_dims[i] << " ";
		}
		os << std::endl;

		// Write the factor scopes
		os << num_factors() << std::endl;
		for (size_t i = 0; i < m_factors.size(); ++i) {
			const factor& f = m_factors[i];
			os << f.nvar();
			variable_set::const_iterator si = f.vars().begin();
			for (; si != f.vars().end(); ++si) {
				os << " " << (*si).label();
			}
			os << std::endl;
		}

		// Write the factor tables
		os << std::endl;
		for (size_t i = 0; i < m_factors.size(); ++i) {
			const factor& f = m_factors[i];
			os << f.numel() << std::endl;
			for (size_t j = 0; j < f.numel(); ++j) {
				os << " " << std::setiosflags(std::ios::fixed)
					<< std::setprecision(8) << f[j];
			}
			os << std::endl << std::endl;
		}
	}

	void write_bayes(std::ostream& os) {

		// Write the header
		size_t n = nvar();
		os << "BAYES" << std::endl;
		os << n << std::endl;
		for (size_t i = 0; i < n; ++i) {
			os << m_dims[i];
			if (i < n - 1) {
				os << " ";
			}
		}
		os << std::endl;

		// Write the factor scopes
		size_t m = m_factors.size();
		std::vector<std::vector<variable> > scopes;
		assert(n == m); // For Bayes nets only
		os << m << std::endl;
		for (size_t i = 0; i < m; ++i) {
			variable_set scope = m_factors[i].vars();
			int child = m_factors[i].get_child();
			std::vector<variable> temp;
			for (variable_set::const_iterator vi = scope.begin();
					vi != scope.end(); ++vi) {
				if (vi->label() != (size_t)child) {
					temp.push_back(*vi);
				}
			}

			temp.push_back(var(child)); // child is last variable
			scopes.push_back(temp);

			os << temp.size();
			for (size_t j = 0; j < temp.size(); ++j) {
				os << " " << temp[j].label();
			}

			os << std::endl;
		}

		os << std::endl;

		// Write the factor tables
		for (size_t i = 0; i < m; ++i) {
			const factor& fn = m_factors[i];
			os << fn.numel() << std::endl;
			std::vector<variable> orig_scope;
			for (variable_set::const_iterator vi = fn.vars().begin();
					vi != fn.vars().end(); ++vi) {
				orig_scope.push_back(*vi);
			}
			convert_index cv(scopes[i], false, orig_scope, true);
			for (size_t j = 0; j < fn.numel(); ++j) {
				size_t k = cv.convert(j);
				os << " " << std::setiosflags(std::ios::fixed)
					<< std::setprecision(8) << fn.get(k);
			}

			os << std::endl << std::endl;
		}

	}
	///
	/// \brief Check if graphical model is a Markov net
	///
	bool is_markov() {
		return m_markov;
	}

	///
	/// \brief Check if the graphical model is a Bayes net
	bool is_bayes() const {
		return !m_markov;
	}

	///
	/// \brief Assert evidence into the graphical model.
	/// \param [in] file_name 		The full path to a file containing the evidence
	/// \param [out] evidence 		The map containing evidence variable value pairs
	/// \param [out] old_to_new 	The mappings between old and new variable indexes
	/// \return the list of new factors resulting from conditioning on the evidence.
	///
	std::vector<factor> assert_evidence(const char* file_name,
			std::map<size_t, size_t>& evidence,
			std::map<size_t, size_t>& old_to_new) {

		std::ifstream is(file_name);
		if (is.fail()) {
			std::cout << "Error while reading evidence file: " << file_name << std::endl;
			throw std::runtime_error("Evidence file error");
		}

		size_t num_vars;
		is >> num_vars; // first line is the number of evidence variable
		for (size_t i = 0; i < num_vars; ++i) {
			size_t var, val;
			is >> var >> val; // pair variable, value
			evidence[var] = val;
		}

		// Condition the factors on the evidence, and accumulate the constant
		variable_set non_evid_vars;
		std::vector<factor> temp;
		std::vector<factor>::const_iterator it = m_factors.begin();
		for (; it != m_factors.end(); ++it) {
			factor f = (*it);
			for (std::map<size_t, size_t>::const_iterator mi = evidence.begin();
					mi != evidence.end(); ++mi) {
				size_t v = mi->first;
				size_t k = mi->second;
				if (f.vars().contains(var(v))) {
					f = f.condition(var(v), k);
				}
			}

			if (f.isscalar()) {
				m_global_const *= f; // accumulate the global constant
			} else {
				temp.push_back(f);
				non_evid_vars |= f.vars();
			}
		}

		// Re-index the non-evidence variables
		vindex idx = 0;
		for (variable_set::const_iterator ci = non_evid_vars.begin();
				ci != non_evid_vars.end(); ++ci) {
			old_to_new[ci->label()] = idx;
			++idx;
		}

		std::vector<factor> result;
		for (size_t i = 0; i < temp.size(); ++i) {
			const factor& f = temp[i];
			variable_set vs;
			for (variable_set::const_iterator ci = f.vars().begin();
					ci != f.vars().end(); ++ci) {
				vindex old_idx = ci->label();
				size_t dim = ci->states();
				vindex new_idx = old_to_new[old_idx];
				vs |= variable(new_idx, dim);
			}

			result.push_back(factor(vs, (double*)f.table()));
		}

		return result;
	}

	///
	/// \brief Assert evidence into the graphical model.
	/// \param [in] evidence 		The map containing evidence variable value pairs
	/// \param [out] old_to_new 	The mappings between old and new variable indexes
	/// \return the list of new factors resulting from conditioning on the evidence.
	///
	std::vector<factor> assert_evidence(const std::map<size_t, size_t>& evidence,
			std::map<size_t, size_t>& old_to_new) {

		// Condition the factors on the evidence, and accumulate the constant
		variable_set non_evid_vars;
		std::vector<factor> temp;
		std::vector<factor>::const_iterator it = m_factors.begin();
		for (; it != m_factors.end(); ++it) {
			factor f = (*it);
			for (std::map<size_t, size_t>::const_iterator mi = evidence.begin();
					mi != evidence.end(); ++mi) {
				size_t v = mi->first;
				size_t k = mi->second;

				// check if value is in domain
				if (k >= var(v).states()) {
					std::stringstream ss;
					ss << "Evidence value " << k << " is outside variable's x"
						<< v << " domain.";
					throw std::runtime_error(ss.str());
				}

				if (f.vars().contains(var(v))) {
					f = f.condition(var(v), k);
				}
			}

			if (f.isscalar()) {
				m_global_const *= f; // accumulate the global constant
			} else {
				temp.push_back(f);
				non_evid_vars |= f.vars();
			}
		}

		// Re-index the non-evidence variables
		vindex idx = 0;
		for (variable_set::const_iterator ci = non_evid_vars.begin();
				ci != non_evid_vars.end(); ++ci) {
			old_to_new[ci->label()] = idx;
			++idx;
		}

		std::vector<factor> result;
		for (size_t i = 0; i < temp.size(); ++i) {
			const factor& f = temp[i];
			variable_set vs;
			for (variable_set::const_iterator ci = f.vars().begin();
					ci != f.vars().end(); ++ci) {
				vindex old_idx = ci->label();
				size_t dim = ci->states();
				vindex new_idx = old_to_new[old_idx];
				vs |= variable(new_idx, dim);
			}

			result.push_back(factor(vs, (double*)f.table()));
		}

		return result;
	}


	// Basic accessors:

	///
	/// \brief Return the number of variables in the model.
	///
	size_t nvar() const {
		return m_vadj.size();
	};

	///
	/// \brief Convert a variable index to Var object.
	/// \param i 	The variable index to be converted
	/// \return the Var object corresponding to the index.
	///
	variable var(vindex i) const {
		return variable(i, m_dims[i]);
	};

	///
	/// \brief Return the number of factors in the model.
	///
	size_t num_factors() const {
		return m_factors.size();
	};

	///
	/// \brief Accessor for a factor.
	/// \param idx 	The index of the factor
	/// \return the factor corresponding to that index.
	///
	const factor& get_factor(findex idx) const {
		return m_factors[idx];
	};

	///
	/// \brief Set a factor
	/// \param i	The index of the factor
	/// \param f	The new factor
	void set_factor(findex i, const factor& f) {
		m_factors[i] = f;
	}

	///
	/// \brief Accessor for the factor container.
	/// \return the list of factors of the model.
	///
	const std::vector<factor>& get_factors() const {
		return m_factors;
	};

	// Basic variable-based queries:

	///
	/// \brief Factors depending on a variable.
	/// \param v 	The variable object
	/// \return the list of factor indexes containing the variable in their scopes.
	///
	const flist& with_variable(const variable& v) const {
		return m_vadj[_vindex(v)];
	}

	///
	/// \brief Union of factors depending on a set of variables.
	/// \param vs 	The set of variable object
	/// \return the list of factor indexes containing those variables in their scopes.
	///
	flist with_variable_set(const variable_set& vs) const { //   or on all of a set of variables
		flist fs = with_variable(vs[0]);
		for (size_t v = 1; v < vs.size(); v++)
			fs &= with_variable(vs[v]);
		return fs;
	}

	
	///
	/// \brief Intersection of factors depending on a set of variables.
	/// \param vs 	The set of variable object
	/// \return the list of factor indexes containing those variables in their scopes.
	///
	flist intersects(const variable_set& vs) const {
		flist fs = with_variable(vs[0]);
		for (size_t v = 1; v < vs.size(); v++)
			fs |= with_variable(vs[v]);
		return fs;
	}

	///
	/// \brief Factors depending on a variable.
	/// \param v 	The variable object
	/// \return the list of factor indexes containing the variable in their scopes.
	///	
	flist contains(const variable& v) const {
		return with_variable(v);
	}

	///
	/// \brief Union of factors depending on a set of variables.
	/// \param vs 	The set of variable object
	/// \return the list of factor indexes containing those variables in their scopes.
	///
	flist contains(const variable_set& vs) const {
		return with_variable_set(vs);
	}
	
	///
	/// \brief Factors contained by a set of variables.
	/// \param vs 	The set of variables
	/// \return the list of factor indexes contained by the set of variables.
	///
	flist contained_by(const variable_set& vs) const {
		flist fs2, fs = intersects(vs);
		for (size_t f = 0; f < fs.size(); f++)
			if (m_factors[fs[f]].vars() << vs)
				fs2 |= fs[f];
		return fs2;
	}

	///
	/// \brief Markov blanket of a variable.
	/// \param v 	The variable object
	/// \return the set of variables that form the Markov blanket 
	/// 	(ie, the variables that *v* may depend on).
	///
	variable_set markov_blanket(const variable& v) const {  // variables that v may depend on
		variable_set vs;
		const flist& nbrs = with_variable(v);
		for (flist::const_iterator f = nbrs.begin(); f != nbrs.end(); ++f)
			vs |= get_factor(*f).vars();
		vs /= v;
		return vs;
	}

	///
	/// \brief Markov blanket of a set of variables.
	/// \param vs 	The set of variables
	/// \return the union of Markov blankets corresponding to the 
	/// 	variables in the input set.
	///
	variable_set markov_blanket(const variable_set& vs) const {
		variable_set ret = markov_blanket(vs[0]);
		for (size_t v = 1; v < vs.size(); v++)
			ret |= markov_blanket(vs[v]);
		return ret;
	}

	///
	/// \brief Full adjacency matrix.
	/// \return the adjacency list associated with each variable in the model.
	///
	std::vector<variable_set> mrf() const {		// full variable-to-var adjacency
		std::vector<variable_set> vvs;
		for (size_t v = 0; v < nvar(); ++v)
			vvs.push_back(markov_blanket(var(v)));
		return vvs;
	}

	// Factor ("node") manipulation operations:

	///
	/// \brief Add a new factor to the model.
	/// \param F 	The factor to be added
	/// \return the index associated with the newly added factor.
	///
	virtual findex add_factor(const factor& F) {         // add a factor to our collection
		const variable_set& v = F.vars();
		findex use = add_node();
		//if (use>=nFactors()) _factors.push_back(F); else _factors[use]=F;
		if (use >= num_factors()) {
			if (m_factors.capacity() > num_factors())
				m_factors.push_back(F);
			else {                           // if we'd need to copy, do it manually
				std::vector<factor> tmp;
				tmp.reserve(2 * num_factors());
				tmp.resize(num_factors() + 1);
				for (size_t i = 0; i < m_factors.size(); ++i)
					tmp[i].swap(m_factors[i]);
				tmp[num_factors()] = F;
				m_factors.swap(tmp);
			}
		} else
			m_factors[use] = F;

		insert(m_vadj, use, v);
		if (m_dims.size() < m_vadj.size())
			m_dims.resize(m_vadj.size(), 0);
		for (variable_set::const_iterator i = v.begin(); i != v.end(); ++i) { // look up dimensions if required
			if (m_dims[_vindex(*i)] == 0)
				m_dims[_vindex(*i)] = i->states();// add if we haven't seen this var
			else if (m_dims[_vindex(*i)] != i->states())	//   or check it against our current states
				throw std::runtime_error(
						"Incompatible state dimension in added factor");
		}
		return use;                                  // return the factor index used
	}

	///
	/// \brief Remove a factor from the model.
	/// \param idx 	The index of the factor to be removed
	/// 
	virtual void remove_factor(findex idx) {        // remove a factor from the collection
		erase(m_vadj, idx, get_factor(idx).vars());		// remove from variable lists
		m_factors[idx] = factor();                       // empty its position
		remove_node(idx);								// and remove the node
	}

	///
	/// \brief Remove all factors.
	///
	void clear_factors() {
		m_factors.clear();
		m_vadj.clear();
		m_vadj.resize(nvar());
		graph::clear();
	};

	///
	/// \brief Find the factor with the smallest scope.
	/// \param fl 	The reference to a list of factors
	/// \return the index of the smallest factor (scope-wise).
	/// 
	findex smallest(const flist& fl) {
		assert(fl.size() > 0);
		findex ret = *fl.begin();
		for (flist::const_iterator f = fl.begin(); f != fl.end(); ++f)
			if (get_factor(ret).nvar() > get_factor(*f).nvar())
				ret = *f;
		return ret;
	}

	///
	/// \brief Find the factor with the largest scope.
	/// \param fl 	The reference to a list of factors
	/// \return the index of the largest factor (scope-wise).
	///
	findex largest(const flist& fl) {
		assert(fl.size() > 0);
		findex ret = *fl.begin();
		for (flist::const_iterator f = fl.begin(); f != fl.end(); ++f)
			if (get_factor(ret).nvar() < get_factor(*f).nvar())
				ret = *f;
		return ret;
	}

	// Check graphical model properties:

	///
	/// \brief Check if a binary model.
	/// \return *true* if all variables have at most two values in 
	/// 	their domains. Otherwise return *false*.
	///
	bool is_binary() const {
		for (size_t i = 0; i < m_dims.size(); ++i) {
			if (m_dims[i] > 2)
				return false;
		}
		return true;
	}
	
	///
	/// \brief Check if a pairwise model.
	/// \return *true* if all factors involve at most two variables.
	/// 	Otherwise return *false*.
	///
	bool is_pairwise() const {
		for (size_t i = 0; i < num_factors(); ++i) {
			if (m_factors[i].nvar() > 2)
				return false;
		}
		return true;
	}

	// Distribution-based operators:

	///
	/// \brief Directly calculate joint distribution function.
	/// \param maxsize 	The maximum size allowed for the joint distribution
	/// \return the factor representing the full joint distribution of the model.
	factor joint(size_t maxsize = 0) const {
		if (maxsize) {
			size_t D = 1;
			for (size_t i = 0; i < nvar(); i++) {   // check for joint being "small"
				D *= (m_dims[i] ? m_dims[i] : 1);
				if (D > maxsize)
					throw std::runtime_error("graphical model::joint too large");
			}
		}
		factor F;				// brute force construction of joint table
		for (size_t i = 0; i < num_factors(); i++)
			F *= m_factors[i];
		return F;
	}

	// Ordering: variable (elimination) orders and factor orders

	///
	/// \brief Find the induced width (complexity) of an elimination order.
	/// \param order 	The variable elimination order
	/// \return the induced width of the elimination order.
	///
	size_t induced_width(const variable_order_t& order) const {
		//return pseudo_tree_size(order).first;

		size_t width = 0;
		std::vector<variable_set> adj = mrf();
		size_t n = order.size();
		graph g(n); // create the graph
		for (size_t i = 0; i < adj.size(); ++i) {
			const variable_set& vi = adj[i];
			for (variable_set::const_iterator cj = vi.begin();
					cj != vi.end(); ++cj) {
				size_t j = _vindex(*cj);
				g.add_edge(i, j);
			}
		}

		std::vector<size_t> position(n);
		for (size_t i = 0; i < n; ++i) {
			position[order[i]] = i;
		}

		// eliminate variables and create induced edges
		for (size_t i = 0; i < n; ++i) {
			size_t var = order[i];
			size_t pos = position[var];

			// find the neighbors appearing later in the ordering
			const set<edge_id> ns = g.neighbors(var);
			std::set<size_t> S;
			for (set<edge_id>::const_iterator si = ns.begin();
					si != ns.end(); ++si) {
				size_t j = si->second;
				if (position[j] > pos) {
					S.insert(j);
				}
			}

			// connect the neighbors appearing later in the ordering
			width = std::max(width, S.size());
			std::set<size_t>::iterator it1, it2;
			for (it1 = S.begin(); it1 != S.end(); ++it1) {
				it2 = it1;
				while (++it2 != S.end()) {
					size_t u = *it1, v = *it2;
					if (u != v) g.add_edge(u, v);
				}
			}
		}

		return width;
	}

	///
	/// \brief Find induced width and pseudo tree height of an elimination order.
	/// \param order The reference of a variable elimination order
	/// \return the pair representing the induced width and height of the 
	/// 	pseudo tree corresponding to the variable elimination order.
	/// 
	std::pair<size_t, size_t> pseudo_tree_size(const variable_order_t& order) const {
		size_t width = 0, maxHeight = 0;
		std::vector<size_t> heights;
		heights.resize(nvar(), 0);
		std::vector<variable_set> adj = mrf();

		// Eliminate in the given order of labels, tracking adjacency
		for (variable_order_t::const_iterator i = order.begin(); i != order.end(); ++i) {
			if (*i < 0 || *i >= nvar())
				continue;
			variable_set vi = adj[_vindex(var(*i))];
			for (variable_set::const_iterator j = vi.begin(); j != vi.end(); ++j) {
				adj[_vindex(*j)] |= vi;
				adj[_vindex(*j)] -= variable_set(var(*i), *j);
				heights[_vindex(*j)] = std::max(heights[_vindex(*j)],
						heights[_vindex(var(*i))] + 1);
				maxHeight = std::max(maxHeight, heights[_vindex(*j)]);
				width = std::max(width, adj[_vindex(*j)].nvar());
			}
		}
		return std::make_pair(width, maxHeight);
	}

	///
	/// \brief Find the pseudo tree corresponding to an elimination order.
	/// \param order 	The reference of a variable elimination order
	/// \return the pseudo tree corresponding to the elimination order
	//		represented by a vector containing the parent of each of 
	//		the variables in the model.
	///
	std::vector<vindex> pseudo_tree(const variable_order_t& order) const {
		std::vector<vindex> parents(nvar(), -1);
		std::vector<variable_set> adj = mrf();

		// Eliminate in the given order of labels, tracking adjacency
		for (variable_order_t::const_iterator i = order.begin(); i != order.end(); ++i) {
			if (*i < 0 || *i >= nvar())
				continue;
			variable_set vi = adj[_vindex(var(*i))];
			for (variable_set::const_iterator j = vi.begin(); j != vi.end(); ++j) {
				adj[_vindex(*j)] |= vi;
				adj[_vindex(*j)] -= variable_set(var(*i), *j);
			}

			for (variable_order_t::const_iterator k = i; k != order.end(); ++k) {
				if (vi.contains(var(*k))) {
					parents[*i] = *k;
					break;
				}
			}
		}
		return parents;
	}

	///
	/// \brief Variable ordering methods.
	///
	MER_ENUM( OrderMethod , MinFill,WtMinFill,MinWidth,WtMinWidth,Random );

    ///
    /// \brief Find a variable elimination order.
    /// \param ord_type 	The ordering method
    /// \return the variable ordering corresponding to the method, such that
    ///		the first variable in the ordering is eliminated first.
    ///
	variable_order_t order(OrderMethod ord_type) const {
		variable_order_t order;
		order.resize(nvar());

		if (ord_type == OrderMethod::Random) {	// random orders are treated here
			for (size_t i = 0; i < nvar(); i++)
				order[i] = var(i).label();	// build a list of all the variables
			std::random_shuffle(order.begin(), order.end());// and randomly permute them
			return order;						    		// then return
		}

		std::vector<variable_set> adj = mrf();
		typedef std::pair<double, size_t> node_score_t;
		typedef std::multimap<double, size_t> score_map_t;
		score_map_t scores;
		std::vector<score_map_t::iterator> reverse(nvar());

		for (size_t v = 0; v < nvar(); v++) 				// get initial scores
			reverse[v] = scores.insert(node_score_t(order_score(adj, v, ord_type), v));

		for (size_t ii = 0; ii < nvar(); ++ii) {// Iterate through, selecting variables
			score_map_t::iterator first = scores.begin();// Choose a random entry from among the smallest
			score_map_t::iterator last = scores.upper_bound(first->first);
			std::advance(first, randi(std::distance(first, last)));
			size_t i = first->second;

			order[ii] = var(i).label();  // save its label in the ordering
			scores.erase(reverse[i]);	 // remove it from our list
			variable_set vi = adj[i]; // go through adjacent variables (copy: adj may change)
			variable_set fix;		  //  and keep track of which need updating
			for (variable_set::const_iterator j = vi.begin(); j != vi.end(); ++j) {
				size_t v = _vindex(*j);
				adj[v] |= vi;             // and update their adjacency structures
				adj[v] /= var(i);
				if (fix.size() < scores.size()) {
					if (ord_type == OrderMethod::MinWidth
							|| ord_type == OrderMethod::WtMinWidth)
						fix |= adj[v]; //var(v);				// (width methods only need v, not nbrs)
					else
						fix |= adj[v];	// come back and recalculate their scores
				}
			}
			for (variable_set::const_iterator j = fix.begin(); j != fix.end(); ++j) {
				size_t jj = j->label();
				scores.erase(reverse[jj]);	// remove and update (score,index) pairs
				reverse[jj] = scores.insert(node_score_t(order_score(adj, jj, ord_type), jj));
			}
		}
		return order;
	}

    ///
    /// \brief Find a variable elimination order.
    /// \param ord_type 	The ordering method
    /// \return the variable ordering corresponding to the method, such that
    ///		the first variable in the ordering is eliminated first.
    ///
	variable_order_t order2(OrderMethod ord_type) const {

		// variable order to be computed
		variable_order_t order;

		if (ord_type == OrderMethod::Random) {	// random orders are treated here
			order.resize(nvar());
			for (size_t i = 0; i < nvar(); i++)
				order[i] = var(i).label();	//   build a list of all the variables
			std::random_shuffle(order.begin(), order.end());//   and randomly permute them
			return order;											//   then return
		}

		// create a temporary graph
		graph G(nvar());
		std::vector<variable_set> adj = mrf();
		for (size_t i = 0; i < nvar(); ++i) {
			variable_set& vs = adj[i];
			for (variable_set::const_iterator vi = vs.begin();
					vi != vs.end(); ++vi) {
				size_t j = vi->label();
				if (j != i) G.add_edge(i, j);
			}
		}

		// initialize the order
		size_t n = nvar();
		order.reserve(n);

		// keeps track of node scores
		std::vector<int> scores(n);
		for (size_t i = 0; i < n; ++i) {
			scores[i] = order_score(G, i, ord_type);
		}

		// eliminate nodes until all gone
		int min_score = -1, num_nodes = n;
		size_t width = 0;
		while (num_nodes != 0) {

			// keeps track of minimal score nodes
			std::vector<size_t> candidates; // minimal score of 1 or higher
			std::vector<size_t> simplicial; // simplicial nodes (score 0)

			min_score = std::numeric_limits<int>::max();

			// find node to eliminate
			for (size_t i = 0; i < n; ++i) {
				if (scores[i] == 0) { // score 0
					simplicial.push_back(i);
				} else if (scores[i] < min_score) { // new, lower score (but greater 0)
					min_score = scores[i];
					candidates.clear();
					candidates.push_back(i);
				} else if (scores[i] == min_score) { // current min. found again
					candidates.push_back(i);
				}
			}

			// eliminate all nodes with score=0 -> no edges will have to be added
			for (std::vector<size_t>::iterator it = simplicial.begin();
					it != simplicial.end(); ++it) {
				size_t v = (*it);
				order.push_back(v);
				--num_nodes;
				set<edge_id> temp = G.neighbors(v);
				width = std::max(width, temp.size());
				for (set<edge_id>::const_iterator ci = temp.begin();
						ci != temp.end(); ++ci) {
					size_t i = ci->first, j = ci->second;
					G.remove_edge(i, j);
				}
				G.remove_node(v); // and adj edges
				scores[v] = std::numeric_limits<int>::max();
			}

			// anything left to eliminate? If not, we are done!
			if (min_score == std::numeric_limits<int>::max()) {
				break;
			}

			// Pick one of the minimal score nodes (with score >= 1),
			// breaking ties randomly
			size_t cand = candidates[randi(candidates.size())];
			//size_t cand = candidates[0]; // first candidate in list (lexicograhically)
			order.push_back(cand);
			--num_nodes;

			// remember it's neighbors, to be used later
			set<edge_id> nlist = G.neighbors(cand);
			std::set<size_t> neighbors;
			for (set<edge_id>::const_iterator ci = nlist.begin();
					ci != nlist.end(); ++ci) {
				if (ci->first == cand && ci->second != cand)
					neighbors.insert(ci->second);
			}

			// connect neighbors in primal graph
			std::set<size_t>::iterator it1, it2;
			for (it1 = neighbors.begin(); it1 != neighbors.end(); ++it1) {
				it2 = it1;
				while (++it2 != neighbors.end()) {
					size_t i = *it1, j = *it2;
					G.add_edge(i, j);
				}
			}

			// compute candidates for score update (node's neighbors and their neighbors)
			width = std::max(width, neighbors.size());
			std::set<size_t> to_fix(neighbors);
			for (std::set<size_t>::const_iterator it = neighbors.begin();
					it != neighbors.end(); ++it) {
				size_t v = (*it);
				const set<edge_id>& temp = G.neighbors(v);
				for (set<edge_id>::const_iterator ci = temp.begin();
						ci != temp.end(); ++ci) {
					if (ci->first == v && v != ci->second)
						to_fix.insert(ci->second);
				}
			}
			to_fix.erase(cand);

			// remove node from primal graph
			for (set<edge_id>::const_iterator ci = nlist.begin();
					ci != nlist.end(); ++ci) {
				size_t i = ci->first, j = ci->second;
				G.remove_edge(i, j);
			}
			G.remove_node(cand); // and adj edges
			scores[cand] = std::numeric_limits<int>::max(); // tag score

			// update scores in primal graph (candidate nodes computed earlier)
			for (std::set<size_t>::const_iterator it = to_fix.begin();
					it != to_fix.end(); ++it) {
				size_t v = (*it);
				scores[v] = order_score(G, v, ord_type);
			}
		}

		return order;
	}

    ///
    /// \brief Find a constrained variable elimination order.
    /// \param ord_type		The ordering method
    /// \param var_types 	The vector containing the variable types (SUM or MAX)
    /// \return the constrained variable ordering corresponding to the method,
    /// such that SUM variables are eliminated before any of the MAX variables.
    ///
	variable_order_t order2(OrderMethod ord_type, std::vector<bool> var_types) const {

		// variable order to be computed
		variable_order_t order;

		// create a temporary graph
		graph G(nvar());
		std::vector<variable_set> adj = mrf();
		for (size_t i = 0; i < nvar(); ++i) {
			variable_set& vs = adj[i];
			for (variable_set::const_iterator vi = vs.begin();
					vi != vs.end(); ++vi) {
				size_t j = _vindex(*vi);
				G.add_edge(i, j);
			}
		}

		// reserve space for the order
		size_t n = nvar();
		order.reserve(n);

		// identify SUM and MAP variables
		std::set<size_t> sum_vars, map_vars;
		for (size_t i = 0; i < n; ++i) {
			if (var_types[i] == false) sum_vars.insert(i);
			else map_vars.insert(i);
		}

		// keeps track of node scores
		std::vector<int> scores(n);
		for (size_t i = 0; i < n; ++i) {
			scores[i] = order_score(G, i, ord_type);
		}

		// eliminate nodes until all gone
		int min_score = -1, num_nodes = n;
		size_t width = 0;
		while (num_nodes != 0) {

			// keeps track of minimal score nodes
			std::vector<size_t> candidates; // minimal score of 1 or higher
			std::vector<size_t> simplicial; // simplicial nodes (score 0)

			min_score = std::numeric_limits<int>::max();

			// check if there are SUM vars left to eliminate
			std::set<size_t>::iterator start, stop;
			if (sum_vars.empty() == false) { // yes
				start = sum_vars.begin();
				stop = sum_vars.end();
			} else { // done, move to MAP vars
				start = map_vars.begin();
				stop = map_vars.end();
			}
			for (std::set<size_t>::iterator si = start; si != stop; ++si) {
				int v = *si;
				if (scores[v] == 0) { // score 0
					simplicial.push_back(v);
				} else if (scores[v] < min_score) { // new, lower score (but greater 0)
					min_score = scores[v];
					candidates.clear();
					candidates.push_back(v);
				} else if (scores[v] == min_score) { // current min. found again
					candidates.push_back(v);
				}
			}

			// eliminate all nodes with score=0 -> no edges will have to be added
			for (std::vector<size_t>::iterator it = simplicial.begin();
					it != simplicial.end(); ++it) {
				size_t v = (*it);
				order.push_back(v);
				--num_nodes;
				set<edge_id> temp = G.neighbors(v);
				width = std::max(width, temp.size());
				for (set<edge_id>::const_iterator ci = temp.begin();
						ci != temp.end(); ++ci) {
					size_t i = ci->first, j = ci->second;
					G.remove_edge(i, j);
				}
				G.remove_node(v); // and adj edges
				scores[v] = std::numeric_limits<int>::max();
				if (var_types[v] == false) sum_vars.erase(v);
				else map_vars.erase(v);
			}

			// anything left to eliminate? If not, we are done!
			if (min_score == std::numeric_limits<int>::max()) {
				break;
			}

			// Pick one of the minimal score nodes (with score >= 1),
			// breaking ties randomly
			size_t cand = candidates[randi2(candidates.size())];
			//size_t cand = candidates[0]; // first candidate in list (lexicograhically)
			order.push_back(cand);
			--num_nodes;

			// remember it's neighbors, to be used later
			set<edge_id> nlist = G.neighbors(cand);
			std::set<size_t> neighbors;
			for (set<edge_id>::const_iterator ci = nlist.begin();
					ci != nlist.end(); ++ci) {
				if (ci->first == cand && ci->second != cand)
					neighbors.insert(ci->second);
			}

			// connect neighbors in primal graph
			std::set<size_t>::iterator it1, it2;
			for (it1 = neighbors.begin(); it1 != neighbors.end(); ++it1) {
				it2 = it1;
				while (++it2 != neighbors.end()) {
					size_t i = *it1, j = *it2;
					G.add_edge(i, j);
				}
			}

			// compute candidates for score update (node's neighbors and their neighbors)
			width = std::max(width, neighbors.size());
			std::set<size_t> to_fix(neighbors);
			for (std::set<size_t>::const_iterator it = neighbors.begin();
					it != neighbors.end(); ++it) {
				size_t v = (*it);
				const set<edge_id>& temp = G.neighbors(v);
				for (set<edge_id>::const_iterator ci = temp.begin();
						ci != temp.end(); ++ci) {
					if (ci->first == v && v != ci->second)
						to_fix.insert(ci->second);
				}
			}
			to_fix.erase(cand);

			// remove node from primal graph
			for (set<edge_id>::const_iterator ci = nlist.begin();
					ci != nlist.end(); ++ci) {
				size_t i = ci->first, j = ci->second;
				G.remove_edge(i, j);
			}
			G.remove_node(cand); // and adj edges
			scores[cand] = std::numeric_limits<int>::max(); // tag score
			if (var_types[cand] == false) sum_vars.erase(cand);
			else map_vars.erase(cand);

			// update scores in primal graph (candidate nodes computed earlier)
			for (std::set<size_t>::const_iterator it = to_fix.begin();
					it != to_fix.end(); ++it) {
				size_t v = (*it);
				scores[v] = order_score(G, v, ord_type);
			}
		}

		return order;
	}

    ///
    /// \brief Find a constrained variable elimination order.
    /// \param ord_type		The ordering method
    /// \param var_types 	The vector containing the variable types (SUM or MAX)
    /// \return the constrained variable ordering corresponding to the method, 
    /// such that SUM variables are eliminated before any of the MAX variables.
    ///
	variable_order_t order(OrderMethod ord_type, std::vector<bool> var_types) const {
		variable_order_t order;
		order.resize(nvar());

		if (ord_type == OrderMethod::Random) {	// random orders are treated here
			variable_order_t sumOrd, maxOrd;
			for (size_t i = 0; i < nvar(); i++) {
				//   build a list of all the variables
				if (var_types[i] == true) maxOrd.push_back(var(i).label());
				else sumOrd.push_back(var(i).label());
			}

			std::random_shuffle(sumOrd.begin(), sumOrd.end());//   and randomly permute them
			std::random_shuffle(maxOrd.begin(), maxOrd.end());
			size_t i = 0;
			for (size_t j = 0; j < sumOrd.size(); ++j) order[i++] = sumOrd[j];
			for (size_t j = 0; j < maxOrd.size(); ++j) order[i++] = maxOrd[j];
			return order;											//   then return
		}

		std::vector<variable_set> adj = mrf();
		typedef std::pair<double, size_t> node_score_t;
		typedef std::multimap<double, size_t> score_map_t;
		score_map_t scores;
		std::vector<score_map_t::iterator> reverse(nvar());

		// get the sum and max variables
		variable_set sumVars;
		for (size_t v = 0; v < nvar(); ++v) {
			if (var_types[v] == false) sumVars |= var(v);
		}

		// get initial scores
		for (size_t v = 0; v < nvar(); ++v) {
			reverse[v] = scores.insert(node_score_t(order_score(adj, v, ord_type), v));
		}

		// Iterate through, selecting sum variables first
		for (size_t ii = 0; ii < nvar(); ++ii) {// Iterate through, selecting variables

			score_map_t::iterator first, last;
			if (sumVars.size() > 0) { // Choose a random entry from among the smallest sum vars
				for (score_map_t::iterator fi = scores.begin(); fi != scores.end(); ++fi) {
					variable v = var(fi->second);
					if (sumVars.contains(v)) {
						first = fi;
						last = scores.upper_bound(first->first);
						break;
					}
				}
			} else {
				first = scores.begin();// Choose a random entry from among the smallest
				last = scores.upper_bound(first->first);
				std::advance(first, randi(std::distance(first, last)));
			}
			size_t i = first->second;

			order[ii] = var(i).label(); // save its label in the ordering
			scores.erase(reverse[i]);	// remove it from our list
			sumVars /= var(i);
			variable_set vi = adj[i]; // go through adjacent variables (copy: adj may change)
			variable_set fix;					//   and keep track of which need updating
			for (variable_set::const_iterator j = vi.begin(); j != vi.end(); ++j) {
				size_t v = _vindex(*j);
				adj[v] |= vi;             // and update their adjacency structures
				adj[v] /= var(i);
				if (fix.size() < scores.size()) {
					if (ord_type == OrderMethod::MinWidth
							|| ord_type == OrderMethod::WtMinWidth)
						fix |= adj[v]; //var(v);				// (width methods only need v, not nbrs)
					else
						fix |= adj[v];	// come back and recalculate their scores
				}
			}
			for (variable_set::const_iterator j = fix.begin(); j != fix.end(); ++j) {
				size_t jj = j->label();
				scores.erase(reverse[jj]);	// remove and update (score,index) pairs
				reverse[jj] = scores.insert(node_score_t(order_score(adj, jj, ord_type), jj));
			}
		}

		return order;
	}

	///
	/// \brief Look for a better order iteratively.
	/// \param ord_type 	The variable ordering method
	/// \param old_score	The best score found so far
	/// \param old_order	The order corresponding to the best score so far
	///
	void improve_order(OrderMethod ord_type, double& old_score, variable_order_t& old_order) const {
		variable_order_t order;
		order.resize(nvar());

		if (ord_type == OrderMethod::Random) {	// random orders are treated here
			for (size_t i = 0; i < nvar(); i++)
				order[i] = var(i).label();	//   build a list of all the variables
			std::random_shuffle(order.begin(), order.end());//   and randomly permute them
			// !!! what scoring mechanism to use?
			//std::pair<size_t,size_t> newsize = pseudoTreeSize(order);
			//if (newsize.first < width || (newsize.first == width && newsize.second < height)) {
			//  width=newsize.first; height=newsize.second; oldOrder=order;
			//}
			old_order = order;
			return;
		}

		std::vector<variable_set> adj = mrf();
		typedef std::pair<double, size_t> NN;
		typedef std::multimap<double, size_t> sMap;
		sMap scores;
		std::vector<sMap::iterator> reverse(nvar());
		double newScore;

		for (size_t v = 0; v < nvar(); v++) 				// get initial scores
			reverse[v] = scores.insert(NN(order_score(adj, v, ord_type), v));

		for (size_t ii = 0; ii < nvar(); ++ii) {// Iterate through, selecting variables
			sMap::iterator first = scores.begin();// Choose a random entry from among the smallest
			newScore = std::max(first->first, newScore);// keep track of largest score for this ordering
			if (newScore > old_score)
				return;            //   & if worse than cutoff value, just quit
			sMap::iterator last = scores.upper_bound(first->first);
			std::advance(first, randi(std::distance(first, last)));
			size_t i = first->second;

			order[ii] = var(i).label();  	       // save its label in the ordering
			scores.erase(reverse[i]);					// remove it from our list
			variable_set vi = adj[i]; // go through adjacent variables (copy: adj may change)
			variable_set fix;					//   and keep track of which need updating
			for (variable_set::const_iterator j = vi.begin(); j != vi.end(); ++j) {
				size_t v = _vindex(*j);
				adj[v] |= vi;             // and update their adjacency structures
				adj[v] /= var(i);
				if (ord_type == OrderMethod::MinWidth
						|| ord_type == OrderMethod::WtMinWidth)
					fix |= adj[v]; //var(v);				// (width methods only need v, not nbrs)
				else
					fix |= adj[v];		// come back and recalculate their scores
			}
			for (variable_set::const_iterator j = fix.begin(); j != fix.end(); ++j) {
				size_t jj = j->label();
				scores.erase(reverse[jj]);	// remove and update (score,index) pairs
				reverse[jj] = scores.insert(NN(order_score(adj, jj, ord_type), jj));
			}
		}
		// !!! check for tiebreaker?
		old_order = order;
		return;
	}

	// Helpful manipulation functions for other data:

	///
	/// \brief Add the scope of a factor to the adjacency list
	/// \param adj 	The adjacency list to be modified
	/// \param idx 	The index of the factor
	/// \param vs 	The scope of the factor
	///
	void insert(std::vector<flist>& adj, findex idx, const variable_set& vs) {
		if (vs.nvar() > 0 && adj.size() <= vs.rbegin()->label()) // if we need to, expand our set of variables
			adj.resize(vs.rbegin()->label() + 1); //   to be large enough to be indexed by label
		for (size_t i = 0; i < vs.nvar(); ++i)
			adj[vs[i]] |= idx;			//   and add factor to adj list
	}

	///
	/// \brief Remove the scope of a factor from the adjacency list
	/// \param adj 	The adjacency list to be modified
	/// \param idx 	The index of the factor
	/// \param vs 	The scope of the factor
	///
	void erase(std::vector<flist>& adj, findex idx, const variable_set& vs) {
		for (size_t i = 0; i < vs.nvar(); ++i)
			adj[vs[i]] /= idx; 		//   remove a factor from each var's adj list
	}

	// Simple optimum selection routines:

	///
	/// \brief Find smallest factor with each variable and select its optimum.
	///
	std::vector<size_t> max_simple() {
		std::vector<size_t> vals(nvar());
		for (size_t v = 0; v < nvar(); ++v)
			if (with_variable(var(v)).size()) {
				vals[v] = get_factor(smallest(with_variable(var(v)))).argmax();
			}
		return vals;
	}

	///
	/// \brief Find *maxmarginal* for variable and select optimum given those so far.
	///
	std::vector<size_t> max_sequential(const variable_order_t& order) {
		std::vector<size_t> vals(nvar());
		variable_set done;
		for (variable_order_t::const_reverse_iterator v = order.rbegin(); v != order.rend();
				++v) {
			variable V = var(*v);
			const flist& use = with_variable(V);
			factor mm(V, 1.0);
			for (flist::const_iterator f = use.begin(); f != use.end(); ++f) {
				size_t condi = 0;
				variable_set condv = get_factor(*f).vars() & done;
				if (condv.size())
					condi = sub2ind(condv, vals);
				mm *= get_factor(*f).condition(condv, condi).maxmarginal(V);
			}
			vals[*v] = mm.argmax();
			done += V;
		}
		return vals;
	}

	///
	/// \brief Compute the log-probability of a particular (complete) configuration.
	///
	double logP(const std::vector<size_t>& config) {
		double val = 0.0;
		for (size_t f = 0; f < num_factors(); ++f)
			val += std::log(get_factor(f)[sub2ind(get_factor(f).vars(), config)]);
		return val;
	}

	///
	/// \brief Output operator.
	/// \param out 	The reference of an output stream
	/// \param gm 	The reference of a graphical model
	/// \return the reference of the modified output stream containing the 
	/// 	formatted content of the graphical model.
	///
	friend std::ostream& operator<<(std::ostream& out, const graphical_model& gm) {
		out << "Dumping graphical model content with " << gm.nvar() << " variables and "
				<< gm.num_factors() << " factors: " << std::endl;
		for (size_t j = 0; j < gm.get_factors().size(); j++)
			out << " " << j << " " << gm.get_factors()[j] << std::endl;
		return out;
	}

	///
	/// \brief Return the constant resulted from asserting evidence in the model.
	///
	double get_global_const() const {
		return m_global_const.max();
	}

	///
	/// \brief Initialize the CPTs uniformly (Bayesian networks only)
	///
	void uniform_bayes() {
		if (is_bayes()) {
			for (size_t i = 0; i < m_factors.size(); ++i) {
				m_factors[i].fill_uniform_bayes();
			}
		}
	}

	///
	/// \brief Initialize the CPTs randomly (Bayesian networks only)
	///
	void random_bayes() {
		if (is_bayes()) {
			for (size_t i = 0; i < m_factors.size(); ++i) {
				m_factors[i].fill_random_bayes();
			}
		}
	}

protected:

	///
	/// \brief Look up the index of a variable.
	/// \param v 	The variable object
	/// \return the index corresponding to the variable in the model.
	size_t _vindex(const variable& v) const {
		return v.label();
	}

	///
	/// \brief Mutable accessor of adjacency.
	/// \param v 	The variable object
	/// \return the list of factors containing the variable in their scope.
	///
	flist& _with_variable(const variable& v) {
		return m_vadj[_vindex(v)];
	}

	///
	/// \brief Internal function for ensuring consistency of the model.
	///
	void fixup() {
		size_t nVar = 0;
		m_vadj.clear();
		m_dims.clear();
		for (std::vector<merlin::factor>::iterator f = m_factors.begin();
				f != m_factors.end(); ++f) {
			if (f->nvar())
				nVar = std::max(nVar, f->vars().rbegin()->label() + 1);
		}

		m_vadj.resize(nVar); // reserve the variable adjaceny lists (factors)
		m_dims.resize(nVar); // make space for variable inclusion mapping
		for (size_t f = 0; f < m_factors.size(); ++f) {	// for each factor,
			findex use = add_node(); // add a node in the graph (factor nodes indexed from 0)
			assert(use == f);
			const variable_set& v = m_factors[f].vars(); // save the variables' dimensions and
			for (variable_set::const_iterator i = v.begin(); i != v.end(); ++i) { // index this factor as including them
				m_dims[_vindex(*i)] = i->states(); // check against current values???
				_with_variable(*i) |= f;
			}
		}
	}

	// Internal helper functions:

	///
	/// \brief Compute the score of a variable for the elimination order.
	/// \param adj 		The adjacency lists
	/// \param i 		The index of the variable
	/// \param kOType 	The variable elimination order
	/// \return the score of the variable.
	///
	double order_score(const std::vector<variable_set>& adj, size_t i, OrderMethod kOType) const {
		double s = 0.0;
		switch (kOType) {
		case OrderMethod::MinFill:
			for (variable_set::const_iterator j = adj[i].begin(); j != adj[i].end(); ++j)
				s += (adj[i] - adj[_vindex(*j)]).size();
			break;
		case OrderMethod::WtMinFill:
			for (variable_set::const_iterator j = adj[i].begin(); j != adj[i].end(); ++j)
				s += (adj[i] - adj[*j]).num_states();
			break;
		case OrderMethod::MinWidth:
			s = adj[i].size();
			break;
		case OrderMethod::WtMinWidth:
			s = adj[i].num_states();
			break;
		default:
			throw std::runtime_error("Unknown elimination ordering type");
			break;
		}
		return s;
	}

	int order_score(const graph& g, size_t v, OrderMethod kOType) const {
		int s = 0.0;
		switch (kOType) {
		case OrderMethod::MinFill:
			{
				const set<edge_id>& nlist = g.neighbors(v);
				std::set<size_t> S;
				for (set<edge_id>::const_iterator ci = nlist.begin(); ci != nlist.end(); ++ci) {
					size_t i = ci->first, j = ci->second;
					if (i == v && i != j)
						S.insert(j);
				}
				std::set<size_t>::iterator it1, it2;
				for (it1 = S.begin(); it1 != S.end(); ++it1) {
					it2 = it1;
					while (++it2 != S.end()) {
						size_t i = *it1, j = *it2;
						if (g.edge(i, j) == edge_id::NO_EDGE)
							++s;
					}
				}
				break;
			}
		case OrderMethod::WtMinFill:
			{
				const set<edge_id>& nlist = g.neighbors(v);
				std::set<size_t> S;
				for (set<edge_id>::const_iterator ci = nlist.begin(); ci != nlist.end(); ++ci) {
					size_t i = ci->first, j = ci->second;
					if (i == v && i != j)
						S.insert(j);
				}
				variable_set vs;
				std::set<size_t>::iterator it1, it2;
				for (it1 = S.begin(); it1 != S.end(); ++it1) {
					it2 = it1;
					while (++it2 != S.end()) {
						size_t i = *it1, j = *it2;
						if (g.edge(i, j) == edge_id::NO_EDGE)
							vs |= var(j);
					}
				}
				s = vs.num_states();
				break;
			}
		case OrderMethod::MinWidth:
			{
				const set<edge_id>& nlist = g.neighbors(v);
				for (set<edge_id>::const_iterator ci = nlist.begin(); ci != nlist.end(); ++ci) {
					size_t i = ci->first, j = ci->second;
					if (i == v && i != j)
						++s;
				}
				break;
			}
		case OrderMethod::WtMinWidth:
			{
				const set<edge_id>& nlist = g.neighbors(v);
				for (set<edge_id>::const_iterator ci = nlist.begin(); ci != nlist.end(); ++ci) {
					size_t i = ci->first, j = ci->second;
					if (i == v && i != j)
						s += var(j).states();
				}
				break;
			}
		default:
			throw std::runtime_error("Unknown elimination ordering type");
			break;
		}
		return s;
	}

	///
	/// \brief Create a random elimination order.
	///
	variable_order_t order_random() const {
		variable_order_t order;
		order.resize(nvar());
		for (size_t i = 0; i < nvar(); i++)
			order[i] = var(i).label();		// build a list of all the variables
		std::random_shuffle(order.begin(), order.end());// and randomly permute them
		return order;
	}

protected:
	// Members:

	std::vector<factor> m_factors;  ///< Collection of all factors in the model.
	bool m_markov;					///< Markov or Bayes network. Default is true (Markov)

private:
	// Members:

	std::vector<flist> m_vadj;		///< Variable adjacency lists (variables to factors)
	std::vector<double> m_dims;		///< Dimensions of variables as stored in graphical model object
	factor m_global_const;			///< Constant produced by removing evidence (default 1.0)
};

} // namespace


#endif /* IBM_MERLIN_GRAPHICAL_MODEL_H_ */
