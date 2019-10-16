/*
 * jglp.cpp
 *
 *  Created on: 12 Sep 2018
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

/// \file jglp.cpp
/// \brief Join Graph Linear Programming (JGLP) algorithm
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#include "jglp.h"

namespace merlin {

void jglp::init() {

	// Prologue
	std::cout << "[JGLP] + i-bound          : " << m_ibound << std::endl;
	std::cout << "[JGLP] + iterations       : " << m_num_iter << std::endl;
	std::cout << "[JGLP] + inference task   : " << "MAP" << std::endl;
	std::cout << "[JGLP] + ordering heur.   : " << m_order_method << std::endl;
	std::cout << "[JGLP] + elimination      : ";

	m_logz = 0.0;
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
	std::cout << "[JGLP] + induced width    : " << wstar << std::endl;
	std::cout << "[JGLP] + exact inference  : " << (m_ibound >= wstar ? "Yes" : "No") << std::endl;
	std::cout << "[JGLP] + ordering time    : " << (timeSystem() - m_start_time) << " seconds" << std::endl;
	if (m_ibound >= wstar)
		m_num_iter = 1; // exact inference requires 1 iteration over the join-tree
}

void jglp::tighten(size_t nIter, double stopTime, double stopObj) {

	std::cout << "[JGLP] Begin iterative cost-shifting over join graph ..." << std::endl;

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
			m_logz += lnmaxf;
			dObj -= lnmaxf;
		}

		// keep track of the tightest upper bound (and configuration)
		if (m_logz < minZ) {
			minZ = m_logz;
			m_best_config = config();
		}

		std::cout << "  logZ: " << std::fixed
			<< std::setw(12) << std::setprecision(MERLIN_PRECISION)
			<< m_logz << " (" << std::scientific
			<< std::setprecision(MERLIN_PRECISION)
			<< std::exp(m_logz) << ") ";
		std::cout << "\td=" << dObj << "\t time="  << std::fixed
			<< std::setprecision(MERLIN_PRECISION) << (timeSystem() - m_start_time)
			<< "\ti=" << iter << std::endl;

	} // done

	// final logZ is the tightest one
	m_logz = minZ;

	double Zdist = std::exp(m_logz / num_factors());
	for (size_t f = 0; f < num_factors(); ++f)
		m_factors[f] *= Zdist;  // !!!! WEIRD; FOR GLOBAL CONSTANT TRANFER

	// Output solution (UAI output format)
	std::cout << "[JGLP] Converged after " << iter << " iterations in "
			<< (timeSystem() - m_start_time) << " seconds" << std::endl;
}

// Write the solution to the output stream
void jglp::write_solution(std::ostream& out, const std::map<size_t, size_t>& evidence,
		const std::map<size_t, size_t>& old2new, const graphical_model& orig,
		const std::set<size_t>& dummies, int output_format) {

	if (output_format == MERLIN_OUTPUT_JSON) {
		out << "{";
		out << " \"algorithm\" : \"jglp\", ";
		out << " \"ibound\" : " << m_ibound << ", ";
		out << " \"iterations\" : " << m_num_iter << ", ";
		out << " \"task\" : \"MAP\", ";
		out << " \"value\" : " << std::fixed
			<< std::setprecision(MERLIN_PRECISION)
			<< (m_logz + std::log(orig.get_global_const())) << ", ";
		out << " \"status\" : \"true\", ";
		out << " \"solution\" : [ ";

		for (vindex i = 0; i < orig.nvar(); ++i) {
			if (dummies.find(i) != dummies.end()) {
				continue; // skip dummy variables from output
			}

			out << "{";
			out << " \"variable\" : " << i << ",";

			try { // evidence variable
				size_t val = evidence.at(i);
				out << " \"value\" : " << val << ",";
			} catch(std::out_of_range& e) { // non-evidence variable
				vindex j = old2new.at(i);
				out << " \"value\" : " << m_best_config[j];
			}
			out << "}";
			if (i != orig.nvar() - 1) {
				out << ", ";
			}
		}
		out << "] ";
		out << "}";
	} else if (output_format == MERLIN_OUTPUT_UAI) {
		out << "MAP" << std::endl;
		out << orig.nvar() - dummies.size();
		for (vindex i = 0; i < orig.nvar(); ++i) {
			if (dummies.find(i) != dummies.end()) {
				continue; // skip dummy variables from output
			}

			try { // evidence variable
				size_t val = evidence.at(i);
				out << " " << val;
			} catch(std::out_of_range& e) { // non-evidence variable
				vindex j = old2new.at(i);
				out << " " << m_best_config[j];
			}
		}
		out << std::endl;
	} else {
		std::string err_msg("Unknown output format.");
		throw std::runtime_error(err_msg);
	}
}

void jglp::run() {

	// Start the timer and store it
	m_start_time = timeSystem();

	// Initialize
	init();

	// Prepare the buckets
	std::vector<factor> fin(m_gmo.get_factors()); // get the original input factors
	std::vector<double> Norm(m_gmo.num_factors(), 0.0); // and normalize them
	for (size_t i = 0; i < m_gmo.num_factors(); ++i) {
		double mx = fin[i].max();
		fin[i] /= mx;
		Norm[i] = std::log(mx);
		m_logz += Norm[i];
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
	std::cout << "[JGLP] Initialize join graph ..." << std::endl;

	for (variable_order_t::const_iterator x = m_order.begin(); x != m_order.end(); ++x) {

		if (m_debug) {
			std::cout << "  Eliminating variable " << *x << std::endl;
		}

		variable VX = var(*x);
		if (*x >= vin.size() || vin[*x].size() == 0)
			continue;  // check that we have some factors over this variable


		// list of factor IDs contained in this bucket
		flist ids = vin[*x];

		// partition into mini-buckets
		partition(VX, fin, vin, Norm, Orig, New, ids);

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

			// Create new cluster alpha over this set of variables;
			// save function parameters also
			findex alpha = findex(-1);//, alpha2 = findex(-1);
			alpha = add_factor(fin[*i]);
			alphas.push_back(alpha);
			m_mini_buckets[*x] |= alpha;

			fin[*i] = elim(fin[*i], VX);

			m_factors[alpha] /= fin[*i];

			// normalize for numerical stability
			double maxf = fin[*i].max();
			fin[*i] /= maxf;
			maxf = std::log(maxf);
			m_logz += maxf;
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
	}
	/// end for: variable elim order

	if (m_debug) {
		std::cout << "Finished creating the join graph." << std::endl;
	}

	// Compute the initial upper bound
	factor F(0.0);
	for (size_t i = 0; i < fin.size(); ++i)
		F += log(fin[i]);
	assert(F.nvar() == 0);
	m_logz += F.max();

	std::cout << "[JGLP] Finished initialization in " << (timeSystem() - m_start_time) << " seconds" << std::endl;
	std::cout << "[JGLP] Initial Upper Bound is "
		<< std::setprecision(MERLIN_PRECISION)
		<< m_logz << " (" << std::scientific
		<< std::setprecision(MERLIN_PRECISION)
		<< std::exp(m_logz)	<< ")" << std::endl;

	// Followed by iterative cost-shifting tighteing via JG propagation
	tighten(m_num_iter);

	// Output the MAP assignment
	std::cout << "[JGLP] Final Upper Bound is "
		<< std::fixed << std::setprecision(MERLIN_PRECISION)
		<< m_logz << " (" << std::scientific
		<< std::setprecision(MERLIN_PRECISION)
		<< std::exp(m_logz)	<< ")" << std::endl;

	m_lb = m_gmo.logP(m_best_config);
	std::cout << "[JGLP] Final Lower Bound is "
		<< std::fixed << std::setprecision(MERLIN_PRECISION)
		<< m_lb << " (" << std::scientific
		<< std::setprecision(MERLIN_PRECISION)
		<< std::exp(m_lb) << ")" << std::endl;

	std::cout << "MAP" << std::endl;
	std::cout << m_best_config.size() << " ";
	std::copy(m_best_config.begin(), m_best_config.end(),
			std::ostream_iterator<index>(std::cout, " "));
	std::cout << std::endl;
}

std::vector<jglp::index> jglp::config() {

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

void jglp::partition(variable VX, std::vector<factor>& fin, std::vector<flist>& vin,
		std::vector<double>& Norm, std::vector<flist>& Orig,
		std::vector<flist>& New, flist& ids) {

	// Select allocation into buckets
	typedef std::pair<double, sPair> _INS;
	std::multimap<double, sPair> scores;
	std::map<sPair, std::multimap<double, sPair>::iterator> reverseScore;

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
			m_logz += mx;
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


} // end namespace


