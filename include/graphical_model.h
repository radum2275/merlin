/*
 * graphical_model.h
 *
 *  Created on: Feb 8, 2013
 *      Author: radu
 */

/// \file graphical_model.h
/// \brief A probabilistic graphical model
/// \author Radu Marinescu


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
			m_factors(), m_vadj(), m_dims() {
	};

	///
	/// \brief Creates a graphical model by copying an existing one.
	/// \param gm 	An object of the same type
	///
	graphical_model(const graphical_model& gm) :
			graph((graph&) gm), m_factors(gm.m_factors), m_vadj(gm.m_vadj), m_dims(
					gm.m_dims) {
	};

	///
	/// \brief Assignment operator (deep copy).
	/// \param gm 	The graphical model to be copied from
	/// \return a reference to the modified object containing the copied content.
	///
	graphical_model& operator=(const graphical_model& gm) {
		graph::operator=((graph&) gm);			// copy graph elements over
		m_factors = gm.m_factors;						// copy list of factors
		m_vadj = gm.m_vadj;								// variable dependencies
		m_dims = gm.m_dims;							// and variable dimensions over
		return *this;
	};

	///
	/// \brief Clone the graphical model.
	/// \return the pointer to the new object representing the copy of 
	/// the current model.
	///
	virtual graphical_model* clone() {
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
	void read(const char* file_name) {

		// Open the input stream
		std::ifstream is(file_name);
		if (is.fail()) {
			std::cout << "Error while opening the input file: " << file_name << std::endl;
			throw std::runtime_error("Input file error");
		}

		// Read the header
		size_t nvar, ncliques, csize, v, nval;
		char st[1024];
		is >> st;
		if ( strcasecmp(st,"MARKOV") )
			throw std::runtime_error("Only UAI Markov-format files are supported currently");

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
		std::vector<factor> tables(ncliques);
		for (size_t i = 0; i < ncliques; i++) {
			is >> nval;
			assert(nval == sets[i].num_states());
			tables[i] = factor(sets[i], 0.0); // preallocate memory and convert from given order, bigEndian
			permute_index pi(cliques[i], true);
			//pi = pi.inverse();   // to our order
			for (size_t j = 0; j < nval; j++) {
				size_t k = pi.convert(j);	// get the index in the factor table
				is >> fval; // read the factor value;
				if (fval == 0.0) fval = 1e-06; // for numerical stability
				tables[i][k] = fval; // save the factor value into the table
			}
		}

		m_factors = tables;
		fixup();
	}

	///
	/// \brief Read a WCNF and convert it to a graphical model.
	/// \param file_name The full path to the input file
	///
	void read_wcnf(const char* file_name) {

		typedef std::pair<double, std::vector<int> > clause_t;

		// Open the input stream
		std::ifstream is(file_name);
		if (is.fail()) {
			std::cout << "Error while opening the input file: " << file_name << std::endl;
			throw std::runtime_error("Input file error");
		}

		// Read the header
		char c;
		std::string wcnf;
		size_t natoms, nclauses; // number of atoms and clauses
		double infinity; // max weight for hard clauses
		is >> c >> wcnf >> natoms >> nclauses >> infinity;
		if ( wcnf.compare("wcnf") != 0 ) {
			throw std::runtime_error("Only WCNF-format files are supported currently");
		}

		// Read the weighted clauses (WCNF: literals are indexed from 1)
		std::vector<bool> unary_hard;
		unary_hard.resize(nclauses, false);
		std::vector<clause_t> clauses;
		for (size_t i = 0; i < nclauses; ++i) {
			double weight; // clause weight
			bool stop = false;
			std::vector<int> lits; // clause literals
			is >> weight; // read the clause weight
			while (!stop) {
				int lit;
				is >> lit; // read the next literal
				if (lit == 0) {
					stop = true; // reached the end of clause (special literal 0)
				} else {
					lits.push_back(lit); // add literal to clause
				}
			}

			// Add the clause to the store
			clauses.push_back(std::make_pair(weight, lits));
			if (weight >= infinity && lits.size() == 1) {
				unary_hard[i] = true;
			}
		}

		// Collect all unique literals and reindex them.
		std::set<int> atoms;
		for (size_t i = 0; i < clauses.size(); ++i) {
			clause_t& cl = clauses[i];
			for (size_t j = 0; j < cl.second.size(); ++j) {
				int lit = cl.second[j];
				atoms.insert(std::abs(lit));
			}
		}

		size_t idx = 0;
		std::map<int, size_t> atom2var;
		for (std::set<int>::iterator si = atoms.begin(); si != atoms.end(); ++si) {
			int lit = (*si);
			atom2var[lit] = idx;
			++idx;
		}

		// Output the atom2var map into a file
		std::string lmap_file_name(file_name);
		lmap_file_name += ".lmap";
		std::ofstream lmap(lmap_file_name.c_str());
		for (std::map<int, size_t>::const_iterator mi = atom2var.begin();
				mi != atom2var.end(); ++mi) {
			lmap << mi->first << " " << mi->second << std::endl;
		}
		lmap.close();

		// Safety checks
		assert(nclauses == clauses.size());
		size_t ncliques = nclauses;
		std::vector<std::vector<variable> > cliques(ncliques);
		std::vector<variable_set> sets(ncliques);

		// Convert each weighted clause into a factor
		std::vector<std::vector<double> > tables(ncliques);
		for (size_t c = 0; c < clauses.size(); ++c) {
			double weight = clauses[c].first;
			std::vector<int> lits = clauses[c].second;
			size_t csize = lits.size();
			cliques[c].reserve(csize);
			std::vector<int> vars;
			vars.resize(csize);
			for (size_t j = 0; j < lits.size(); ++j) {
				int a = std::abs(lits[j]);
				std::map<int, size_t>::iterator mi = atom2var.find(a);
				assert(mi != atom2var.end());
				int v = (*mi).second; // corresp. variable index (from 0)
				variable V(v, 2); // boolean variable
				cliques[c].push_back(V);
				sets[c] |= V;
				vars[j] = v;
			}

			// Loop over the truth assignments of the variables (LSB)
			size_t tsize = sets[c].num_states();

			// Value assignments to the variables
			std::vector<int> vals;
			vals.resize(vars.size(), 0);

			// Reset the assignments to the variables [0, 0, ..., -1]
			vals[csize - 1] = -1;

			// Iterate over the value assignments (least significant digit changes first)
			int i;
			while (true) {

				// find next configuration
				for (i = csize - 1; i >= 0; --i) {
					int var = vars[i];
					int lastVal = 1;
					if (vals[i] < lastVal) break;
					vals[i] = 0;
				}

				if (i < 0) break;	// done;
				++vals[i];

				// Now all variables have a specific value assignment
				// so, compute the corresponding truth assignment of the clause
				bool truth = false;
				for (int j = 0; j < csize; ++j) {
					bool val = (vals[j] == 0 ? false : true);
					if (lits[j] < 0) { // negative literal
						truth |= (!val);
					} else { // positive literal
						truth |= val;
					}

					if (truth) break; // satisfied
				}

				double pval = 1.0; // default false: f(vals) = 1
				if (truth) { // true: f(vals) = e^weight
					if (weight == infinity) pval = infinity;
					else pval = std::exp(weight);
				}
				tables[c].push_back(pval);
			}
		}

		// Create the factors (ensure conversion to ordered scopes)
		std::vector<factor> factors(ncliques);
		for (size_t c = 0; c < ncliques; c++) {
			std::vector<double>& temp = tables[c]; // temporary table
			size_t nval = temp.size();
			assert(nval == sets[c].num_states());
			factors[c] = factor(sets[c], 0.0); // preallocate memory and convert from given order, bigEndian
			permute_index pi(cliques[c], false);
			for (size_t j = 0; j < nval; j++) {
				factors[c][pi.convert(j)] = temp[j];
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
	void write(const char* file_name) {

		// Open the output stream
		std::ofstream os(file_name);
		if (os.fail()) {
			std::cout << "Error while opening the output file: " << file_name << std::endl;
			throw std::runtime_error("Output file error");
		}

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

		// Close the output stream
		os.close();
	}

	///
	/// \brief Write the graphical model to a file using the libDAI .fg format.
	///
	/// For details on the .fg file format see the libDAI documentation. The
	/// internal representation of the factors assume that the scopes are ordered
	/// lexicographically, and the factor table is represented using the "least
	/// significant digit" convention, namely the last variable in the scope
	/// changes its value the fastest (i.e., the UAI convention).
	///
	void write2fg(const char* file_name) {

		// Open the output stream
		std::ofstream os(file_name);
		if (os.fail()) {
			std::cout << "Error while opening the output file: " << file_name << std::endl;
			throw std::runtime_error("Output file error");
		}

		// Write the header and number of factors
		os << "# Factor graph produced by merlin (least significant digit convention)" << std::endl;
		os << num_factors() << std::endl;
		os << std::endl;

		// Write the factors
		for (size_t i = 0; i < m_factors.size(); ++i) {
			const factor& f = m_factors[i];
			os << f.nvar() << std::endl; // scope size
			variable_set::const_iterator si = f.vars().begin(); // scope
			for (; si != f.vars().end(); ++si) {
				os << (*si).label() << " ";
			}
			os << std::endl;
			for (si = f.vars().begin(); si != f.vars().end(); ++si) { // domains
				os << (*si).states() << " ";
			}
			os << std::endl;
			os << f.numel() << std::endl; // table size
			for (size_t j = 0; j < f.numel(); ++j) {
				os << j << " " << std::setiosflags(std::ios::fixed)
					<< std::setprecision(8) << f[j] << std::endl;
			}
			os << std::endl << std::endl;
		}

		// Close the output stream
		os.close();
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

		size_t nvar;
		is >> nvar; // first line is the number of evidence variable
		for (size_t i = 0; i < nvar; ++i) {
			size_t var, val;
			is >> var >> val; // pair variable, value
			evidence[var] = val;
		}

		// Condition the factors on the evidence, and accumulate the constant
		variable_set nonEvidVars;
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
				nonEvidVars |= f.vars();
			}
		}

		// Re-index the non-evidence variables
		vindex idx = 0;
		for (variable_set::const_iterator ci = nonEvidVars.begin();
				ci != nonEvidVars.end(); ++ci) {
			old_to_new[ci->label()] = idx;
			++idx;
		}

		std::vector<factor> result;
		for (size_t i = 0; i < temp.size(); ++i) {
			const factor& f = temp[i];
			variable_set vs;
			for (variable_set::const_iterator ci = f.vars().begin();
					ci != f.vars().end(); ++ci) {
				vindex oldidx = ci->label();
				size_t dim = ci->states();
				vindex newidx = old_to_new[oldidx];
				vs |= variable(newidx, dim);
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
		variable_set nonEvidVars;
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
				nonEvidVars |= f.vars();
			}
		}

		// Re-index the non-evidence variables
		vindex idx = 0;
		for (variable_set::const_iterator ci = nonEvidVars.begin();
				ci != nonEvidVars.end(); ++ci) {
			old_to_new[ci->label()] = idx;
			++idx;
		}

		std::vector<factor> result;
		for (size_t i = 0; i < temp.size(); ++i) {
			const factor& f = temp[i];
			variable_set vs;
			for (variable_set::const_iterator ci = f.vars().begin();
					ci != f.vars().end(); ++ci) {
				vindex oldidx = ci->label();
				size_t dim = ci->states();
				vindex newidx = old_to_new[oldidx];
				vs |= variable(newidx, dim);
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
	flist with_var_set(const variable_set& vs) const { //   or on all of a set of variables
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
		return with_var_set(vs);
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
	findex add_factor(const factor& F) {         // add a factor to our collection
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
	void remove_factor(findex idx) {        // remove a factor from the collection
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
	bool isbinary() const {
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
	bool ispairwise() const {
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
					throw std::runtime_error("graphModel::joint too large");
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
		return pseudo_tree_size(order).first;
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
				order[i] = var(i).label();	//   build a list of all the variables
			std::random_shuffle(order.begin(), order.end());//   and randomly permute them
			return order;											//   then return
		}

		std::vector<variable_set> adj = mrf();
		//vector<set<Var> > adj(adj1.size());
		//for (size_t i=0;i<adj1.size();++i) adj[i] = set<Var>(adj1[i].begin(),adj1[i].end());

		typedef std::pair<double, size_t> NN;
		typedef std::multimap<double, size_t> sMap;
		sMap scores;
		std::vector<sMap::iterator> reverse(nvar());

		for (size_t v = 0; v < nvar(); v++) 				// get initial scores
			reverse[v] = scores.insert(NN(order_score(adj, v, ord_type), v));

		for (size_t ii = 0; ii < nvar(); ++ii) {// Iterate through, selecting variables
			sMap::iterator first = scores.begin();// Choose a random entry from among the smallest
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
				reverse[jj] = scores.insert(NN(order_score(adj, jj, ord_type), jj));
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
		typedef std::pair<double, size_t> NN;
		typedef std::multimap<double, size_t> sMap;
		sMap scores;
		std::vector<sMap::iterator> reverse(nvar());

		// get the sum and max variables
		variable_set sumVars;
		for (size_t v = 0; v < nvar(); ++v) {
			if (var_types[v] == false) sumVars |= var(v);
		}

		// get initial scores
		for (size_t v = 0; v < nvar(); ++v) {
			reverse[v] = scores.insert(NN(order_score(adj, v, ord_type), v));
		}

		// Iterate through, selecting sum variables first
		for (size_t ii = 0; ii < nvar(); ++ii) {// Iterate through, selecting variables

			sMap::iterator first, last;
			if (sumVars.size() > 0) { // Choose a random entry from among the smallest sum vars
				for (sMap::iterator fi = scores.begin(); fi != scores.end(); ++fi) {
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

			order[ii] = var(i).label();  	       // save its label in the ordering
			scores.erase(reverse[i]);					// remove it from our list
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
				reverse[jj] = scores.insert(NN(order_score(adj, jj, ord_type), jj));
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
	flist& _withVariable(const variable& v) {
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
			findex use = add_node();
			assert(use == f);
			const variable_set& v = m_factors[f].vars(); // save the variables' dimensions and
			for (variable_set::const_iterator i = v.begin(); i != v.end(); ++i) { // index this factor as including them
				m_dims[_vindex(*i)] = i->states(); // check against current values???
				_withVariable(*i) |= f;
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

private:
	// Members:

	std::vector<flist> m_vadj;		///< Variable adjacency lists (variables to factors)
	std::vector<double> m_dims;		///< Dimensions of variables as stored in graphical model object
	factor m_global_const;			///< Constant produced by removing evidence (default 1.0)
};

} // namespace

// ToDo:
// Container aspect: contains
//    vector<Factor> : insert(F), erase(i), compress (maximal sets only)
//    map<Var,set<fidx>> : find factors with var / varset
//    edges between factor "nodes" : add, remove, clear, neighborhood lists

// Functional operations
//    Orderings: variable (elimination orders) and node (spanning & covering trees, etc)
//    MRF adj:  map<Var, VarSet> , markovBlanket(v), cliques?
//    VarDepUpdates:  change Var:FIdx map:  erase(&map,i,vs), insert(&map,i,vs)
//    Edge/adjacency:  edge(i,j), neighbors(i),
//




#endif /* GRAPHICAL_MODEL_H_ */
