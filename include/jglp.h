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
/// \author Radu Marinescu


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
	virtual jglp* clone() const {
		jglp* gm = new jglp(*this);
		return gm;
	}

	double ub() const {
		return m_log_z;
	}
	double lb() const {
		return m_lb;
	}
	std::vector<index> best_config() const {
		return m_best_config;
	}

	double logZ() const {
		return m_log_z;
	}
	double logZub() const {
		return m_log_z;
	}
	double logZlb() const {
		return m_log_z;
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
	MER_ENUM( Property , iBound,Order,Iter );

protected:
	// Members:

	graphical_model m_gmo; 					///< Original graphical model
	OrderMethod m_order_method;				///< Variable ordering heuristic
	size_t m_ibound;						///< Mini-bucket i-bound
	double m_log_z;							///< Log partition function
	double m_lb;							///< Lower bound
	variable_order_t m_order; 				///< Elimination ordering
	std::vector<vindex> m_parents; 			///< Pseudo tree
	std::vector<index> m_best_config;		///< MAP assignment
	size_t m_num_iter;						///< Number of iterations
	std::vector<flist> m_mini_buckets;		///< Mini-buucket partitionings

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
			set_properties("iBound=4,Order=MinFill,Iter=100");
			return;
		}
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
	void init() {

		// Start the timer and store it
		m_start_time = timeSystem();

		// Prologue
		std::cout << VERSIONINFO << std::endl << COPYRIGHT << std::endl;
		std::cout << "Initialize inference engine ..." << std::endl;
		std::cout << "+ tasks supported  : MAP" << std::endl;
		std::cout << "+ algorithm        : " << "JGLP" << std::endl;
		std::cout << "+ i-bound          : " << m_ibound << std::endl;
		std::cout << "+ iterations       : " << m_num_iter << std::endl;
		std::cout << "+ inference task   : " << "MAP" << std::endl;
		std::cout << "+ ordering heur.   : " << m_order_method << std::endl;
		std::cout << "+ elimination      : ";

		m_log_z = 0.0;
		if (m_order.size() == 0) { // if we need to construct an elimination ordering
			m_order = m_gmo.order(m_order_method);
			m_parents.clear(); // (new elim order => need new pseudotree) !!! should do together
			std::copy(m_order.begin(), m_order.end(),
					std::ostream_iterator<size_t>(std::cout, " "));
		}
		if (m_parents.size() == 0) {     // if we need to construct a pseudo-tree
			m_parents = m_gmo.pseudo_tree(m_order);
		}

		std::cout << std::endl;
		size_t wstar = m_gmo.induced_width(m_order);
		std::cout << "+ induced width    : " << wstar << std::endl;
		std::cout << "+ exact inference  : " << (m_ibound >= wstar ? "Yes" : "No") << std::endl;
		if (m_ibound >= wstar)
			m_num_iter = 1; // exact inference requires 1 iteration over the join-tree
	}

	///
	/// \brief Join graph propagation for tightening the bound.
	/// \param nIter 	The number of iterations
	/// \param stopTime The time limit
	/// \param stopObj 	The difference between two successive objective values
	///
	void tighten(size_t nIter, double stopTime = -1, double stopObj = -1) {

		std::cout << "Begin iterative cost-shifting over join graph ..." << std::endl;
		std::cout << "- stopObj  : " << stopObj << std::endl;
		std::cout << "- stopTime : " << stopTime << std::endl;
		std::cout << "- stopIter : " << nIter << std::endl;

		double minZ = infty();
		const std::vector<edge_id>& elist = edges();
		double startTime = timeSystem(), dObj = infty();
		size_t iter;
		for (iter = 0; iter < nIter; ++iter) {
			if (std::abs(dObj) < stopObj) {
				break;
			} else {
				dObj = 0.0;
			}

			// iterate over the join-graph edges
			for (size_t i = 0; i < elist.size(); ++i) {
				if (stopTime > 0 && stopTime <= (timeSystem() - startTime)) {
					iter = nIter;
					break;
				}
				findex a, b;
				a = elist[i].first;
				b = elist[i].second;
				if (a > b)
					continue;

				variable_set both = m_factors[a].vars() & m_factors[b].vars();
				//std::cout<<_factors[a].vars()<<"; "<<_factors[b].vars()<<"= "<<both<<"\n";
				factor fratio = (m_factors[a].maxmarginal(both)
						/ m_factors[b].maxmarginal(both)) ^ (0.5);
				m_factors[b] *= fratio;
				m_factors[a] /= fratio;
			}

			// normalize the factors for numerical stability
			for (size_t i = 0; i < num_factors(); ++i) {
				double maxf = m_factors[i].max();
				m_factors[i] /= maxf;
				double lnmaxf = std::log(maxf);
				m_log_z += lnmaxf;
				dObj -= lnmaxf;
			}

			// keep track of the tightest upper bound (and configuration)
			if (m_log_z < minZ) {
				minZ = m_log_z;
				m_best_config = config();
			}

			std::cout << "  JGLP: " << std::fixed << std::setw(12) << std::setprecision(6)
				<< m_log_z << " (" << std::scientific << std::setprecision(6)
				<< std::exp(m_log_z) << ") ";
			std::cout << "\td=" << dObj << "\t time="  << std::fixed
				<< std::setprecision(6) << (timeSystem() - m_start_time)
				<< "\ti=" << iter << std::endl;

		} // done

		// final logZ is the tightest one
		m_log_z = minZ;

		double Zdist = std::exp(m_log_z / num_factors());
		for (size_t f = 0; f < num_factors(); ++f)
			m_factors[f] *= Zdist;  // !!!! WEIRD; FOR GLOBAL CONSTANT TRANFER

		// Output solution (UAI output format)
		std::cout << "Converged after " << iter << " iterations in "
				<< (timeSystem() - m_start_time) << " seconds" << std::endl;
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
			throw std::runtime_error("Error while opening the output file.");
		}

		out << "MAP" << std::endl;
		out << orig.nvar();
		for (vindex i = 0; i < orig.nvar(); ++i) {
			try { // evidence variable
				size_t val = evidence.at(i);
				out << " " << val;
			} catch(std::out_of_range& e) { // non-evidence variable
				vindex j = old2new.at(i);
				out << " " << m_best_config[j];
			}
		}
		out << std::endl;
	}

	void run() {

		// Initialize
		init();

		// Prepare the buckets
		std::vector<factor> fin(m_gmo.get_factors()); // get the original input factors
		std::vector<double> Norm(m_gmo.num_factors(), 0.0); // and normalize them
		for (size_t i = 0; i < m_gmo.num_factors(); ++i) {
			double mx = fin[i].max();
			fin[i] /= mx;
			Norm[i] = std::log(mx);
			m_log_z += Norm[i];
		}

		// Map each variable to the list of factors mentioning that variable
		std::vector<flist> vin;
		for (size_t i = 0; i < m_gmo.nvar(); ++i) {
			vin.push_back(m_gmo.with_variable(var(i)));
		}

		// Origination info: which original factors are included for the first
		// time, and which newly created clusters feed into this cluster
		std::vector<flist> Orig(m_gmo.num_factors());
		std::vector<flist> New(m_gmo.num_factors());
		for (size_t i = 0; i < Orig.size(); ++i) {
			Orig[i] |= i;
		}

		// save the mini-buckets (as lists of factor indeces)
		m_mini_buckets.resize(m_gmo.nvar());

		// Eliminate each variable in the sequence given:
#ifdef DEBUG
		std::cout << "[DEBUG]: Start eliminating variables\n";
#endif

		for (variable_order_t::const_iterator x = m_order.begin(); x != m_order.end(); ++x) {

#ifdef DEBUG
			std::cout << "[DEBUG] Eliminating variable " << *x << " and logZ = " << m_log_z;
#endif

			variable VX = var(*x);
			if (*x >= vin.size() || vin[*x].size() == 0)
				continue;  // check that we have some factors over this variable


			// list of factor IDs contained in this bucket
			flist ids = vin[*x];

			// partition into mini-buckets
			partition(VX, fin, vin, Norm, Orig, New, ids);

#ifdef DEBUG
			std::cout << " with " << ids.size() << " mini-buckets" << "\n";
#endif
			// Perform any matching?
			//    "Matching" here is: compute the largest overlap of all buckets,
			//    and ensure that the
			//    moments on that subset of variables are identical in all buckets.
			if ( ids.size() > 1 ) {
				std::vector<factor> ftmp(ids.size());  // compute geometric mean
				variable_set var = fin[ids[0]].vars();      // on all mutual variables
				for (size_t i = 1; i < ids.size(); i++)
					var &= fin[ids[i]].vars();
				//Factor fmatch(var,1.0);
				factor fmatch(var, 0.0);
				for (size_t i = 0; i < ids.size(); i++) {
					//ftmp[i] = marg(fin[ids[i]],var);
					//fmatch *= ftmp[i];
					ftmp[i] = marg(fin[ids[i]], var).log();
					fmatch += ftmp[i];
				}
				//fmatch ^= (1.0/ids.size());                  // and match each bucket to it
				fmatch *= (1.0 / ids.size());     // and match each bucket to it
				//for (size_t i=0;i<ids.size();i++) fin[ids[i]] *= fmatch/ftmp[i];
				for (size_t i = 0; i < ids.size(); i++)
					fin[ids[i]] *= (fmatch - ftmp[i]).exp();
			}

			// Eliminate individually within buckets
			std::vector<findex> alphas;
			for (flistIt i = ids.begin(); i != ids.end(); ++i) {

#ifdef DEBUG
				std::cout << "[DEBUG] -  processing mb " << *i << " " << fin[*i] << " yields ";
#endif
				// Create new cluster alpha over this set of variables;
				// save function parameters also
				findex alpha = findex(-1);//, alpha2 = findex(-1);
				alpha = add_factor(fin[*i]);
				alphas.push_back(alpha);
				m_mini_buckets[*x] |= alpha;

				fin[*i] = elim(fin[*i], VX);

#ifdef DEBUG
				std::cout << fin[*i] << "\n";
#endif

				m_factors[alpha] /= fin[*i];

				// normalize for numerical stability
				double maxf = fin[*i].max();
				fin[*i] /= maxf;
				maxf = std::log(maxf);
				m_log_z += maxf;
				Norm[*i] += maxf; // save normalization for overall bound

				for (size_t j = 0; j < alphas.size() - 1; ++j)
					add_edge(alpha, alphas[j]);
				for (flistIt j = New[*i].begin(); j != New[*i].end(); ++j)
					add_edge(*j, alpha);

				Orig[*i].clear();
				New[*i].clear();
				New[*i] |= alpha;  // now incoming nodes to *i is just alpha

				insert(vin, *i, fin[*i].vars()); // recompute and update adjacency
			}
			//std::cout<<"\n";

		}
		/// end for: variable elim order

#ifdef DEBUG
		std::cout << "[DEBUG]: Finished eliminating variables\n";
#endif

		factor F(0.0);
		for (size_t i = 0; i < fin.size(); ++i)
			F += log(fin[i]);
		assert(F.nvar() == 0);
		m_log_z += F.max();

		std::cout << "Initial Upper Bound is " << std::fixed << std::setw(12) << std::setprecision(6)
			<< m_log_z << " (" << std::scientific << std::setprecision(6)
			<< std::exp(m_log_z)	<< ")" << std::endl;

		// Followed by iterative cost-shifting tighteing via JG propagation
		tighten(m_num_iter);

		// output the assignment
		std::cout << "Final Upper Bound is " << std::fixed << std::setw(12) << std::setprecision(6)
			<< m_log_z << " (" << std::scientific << std::setprecision(6)
			<< std::exp(m_log_z)	<< ")" << std::endl;

		m_lb = m_gmo.logP(m_best_config);
		std::cout << "Final Lower Bound is " << std::fixed << std::setw(12) << std::setprecision(6)
			<< m_lb << " (" << std::scientific << std::setprecision(6)
			<< std::exp(m_lb) << ")" << std::endl;

		std::cout << "MAP" << std::endl;
		std::cout << m_best_config.size() << " ";
		std::copy(m_best_config.begin(), m_best_config.end(),
				std::ostream_iterator<index>(std::cout, " "));
		std::cout << std::endl;
	}


	///
	/// \brief Return the MAP configuration following the propagation.
	///
	std::vector<index> config() {

		std::vector<index> best;

		// recover the MAP assignment
		best.resize(m_gmo.nvar(), -1);
		variable_set vars;
		for (variable_order_t::reverse_iterator x = m_order.rbegin();
				x != m_order.rend(); ++x) {

			variable VX = var(*x);
			flist ids = m_mini_buckets[*x]; // bucket of VX
			factor F(1.0);
			for (flist::const_iterator i = ids.begin(); i != ids.end(); ++i) {
				factor f = m_factors[*i];
				for (variable_set::const_iterator v = vars.begin(); v != vars.end(); ++v) {
					if (f.vars().contains(*v))
						f = f.condition(variable_set(*v), best[v->label()]);
				}
				F *= f;
			}

			best[*x] = F.argmax();
			vars |= VX;
		}

		return best;
	}

private:

	// Partition the bucket into mini-buckets; initially 'ids' is the bucket,
	// but at the end it contains a list of mini-buckets (factors in each
	// mini-bucket are combined into a single one, thus a mini-bucket is a factor
	// defined over at most iBound distinct variables).
	void partition(variable VX, std::vector<factor>& fin, std::vector<flist>& vin,
			std::vector<double>& Norm, std::vector<flist>& Orig,
			std::vector<flist>& New, flist& ids) {

		// Select allocation into buckets
		typedef std::pair<double, sPair> _INS;
		std::multimap<double, sPair> scores;
		std::map<sPair, std::multimap<double, sPair>::iterator> reverseScore;

//		std::cout<<"Initial table sizes: ";
//			for (flistIt i=ids.begin();i!=ids.end();++i)
//				std::cout<<fin[*i].numel()<<" "; std::cout<<"\n";
//		std::cout << "Initial factors: ";
//		for (flistIt i=ids.begin();i!=ids.end();++i)
//			std::cout<<*i << ": " <<fin[*i]<<"; "; std::cout<<"\n";

		// Populate list of pairwise scores for aggregation
		for (flistIt i = ids.begin(); i != ids.end(); ++i) {
			for (flistIt j = ids.begin(); j != i; ++j) {
				double err = score(fin, VX, *i, *j);
				sPair sp(*i, *j);
				reverseScore[sp] = scores.insert(_INS(err, sp)); // save score
			}
			reverseScore[sPair(*i, *i)] = scores.insert(
					_INS(-1, sPair(*i, *i)));    // mark self index at -1
		}

		// Run through until no more pairs can be aggregated:
		//   Find the best pair (ii,jj) according to the scoring heuristic
		//   and join them as jj; then remove ii and re-score all pairs with jj
		for (;;) {
			std::multimap<double, sPair>::reverse_iterator top = scores.rbegin();
			if (top->first < 0)
				break;                         // if can't do any more, quit
			else {
				size_t ii = top->second.first, jj = top->second.second;
//				std::cout<<"Joining "<<ii<<","<<jj<<"; size "<<(fin[ii].vars()+fin[jj].vars()).nrStates()<<"\n";
				fin[jj] *= fin[ii];                        // combine into j
				Norm[jj] += Norm[ii];
				double mx = fin[jj].max();
				fin[jj] /= mx;
				mx = std::log(mx);
				m_log_z += mx;
				Norm[jj] += mx;
				erase(vin, ii, fin[ii].vars());
				fin[ii] = factor();  //   & remove i

				Orig[jj] |= Orig[ii];
				Orig[ii].clear(); // keep track of list of original factors in this cluster
				New[jj] |= New[ii];
				New[ii].clear(); //   list of new "message" clusters incoming to this cluster

				for (flistIt k = ids.begin(); k != ids.end(); ++k) { // removing entry i => remove (i,k) for all k
					scores.erase(reverseScore[sPair(ii, *k)]);
				}

				ids /= ii;

				for (flistIt k = ids.begin(); k != ids.end(); ++k) { // updated j; rescore all pairs (j,k)
					if (*k == jj)
						continue;
					double err = score(fin, VX, jj, *k);
					sPair sp(jj, *k);
					scores.erase(reverseScore[sp]);    // change score (i,j)
					reverseScore[sp] = scores.insert(_INS(err, sp));  //
				}
			}
		}
	}

};

} // namespace

#endif // re-include
