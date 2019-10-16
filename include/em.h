/*
 * em.h
 *
 *  Created on: Aug 2, 2013
 *      Author: radu
 *
 * Copyright (c) 2019, International Business Machines Corporation
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

/// \file em.h
/// \brief EM parameter learning for Bayes nets (directed models)
/// \author Radu Marinescu radu.marinescu@ie.ibm.com
///

#ifndef IBM_MERLIN_CORE_EM_H_
#define IBM_MERLIN_CORE_EM_H_

#include "base.h"
#include "util.h"
#include "graphical_model.h"
#include "cte.h"

namespace merlin {

class em {

	typedef graphical_model::vindex vindex;
	typedef graphical_model::OrderMethod OrderMethod;
	typedef std::vector<double> likelihood;

	///
	/// \brief Properties of the algorithm.
	///
	MER_ENUM( Property , Order,Infer,Iter,Debug,Threshold,Init );

	///
	/// \brief Inference algorithm.
	///
	MER_ENUM( InferMethod , CTE,WMB );

	///
	/// \brief Factor initialization method.
	///
	MER_ENUM( InitMethod, None,Uniform,Random);

protected:

	///
	/// \brief The expectation (E) step
	///
	void e_step();

	///
	/// \brief The maximization (M) step
	///
	void m_step();

	///
	/// \brief Initialize the EM algorithm
	///
	void init();

public:

	///
	/// \brief Set the training dataset
	/// \param d The dataset
	///
	void set_dataset(std::vector<std::vector<observation> >& d) {
		m_dataset = d;
	}

	///
	/// \brief Set the original graphical model
	/// \param gm The graphical model
	///
	void set_model(graphical_model& gm) {
		m_gmo = gm;
	}

	///
	/// \brief Set the properties of the algorithm.
	/// \param opt 	The string containing comma separated property value pairs
	///
	virtual void set_properties(std::string opt = std::string()) {
		if (opt.length() == 0) {
			set_properties("Order=MinFill,Infer=CTE,Iter=10,Debug=0,Threshold=1e-6,Init=Uniform");
			return;
		}
		m_debug = false;
		std::vector<std::string> strs = merlin::split(opt, ',');
		for (size_t i = 0; i < strs.size(); ++i) {
			std::vector<std::string> asgn = merlin::split(strs[i], '=');
			switch (Property(asgn[0].c_str())) {
			case Property::Order:
				m_order_method = OrderMethod(asgn[1].c_str());
				break;
			case Property::Iter:
				m_iterations = atol(asgn[1].c_str());
				break;
			case Property::Infer:
				m_infer_method = InferMethod(asgn[1].c_str());
				break;
			case Property::Debug:
				m_debug = (atol(asgn[1].c_str()) == 0) ? false : true;
				break;
			case Property::Init:
				m_init_method = InitMethod(asgn[1].c_str());
				break;
			case Property::Threshold:
				m_epsilon = atof(asgn[1].c_str());
				break;
			default:
				break;
			}
		}
	}

	///
	/// \brief Write the solution (i.e., new learned parameters via EM)
	/// \param out The output stream
	/// \param orig	The original graphical model
	///
	void write_solution(std::ostream& out, const graphical_model& orig);

	///
	/// \brief Run the EM algorithm for parameter learning
	///
	void run();

public:

	// Default constructor
	em(graphical_model& gm) : m_iterations(10),
		m_loglikelihood(0),
		m_epsilon(1e-6),
		m_debug(false),
		m_order_method(OrderMethod::MinFill),
		m_infer_method(InferMethod::CTE),
		m_init_method(InitMethod::Uniform),
		m_gmo(gm) {};

	// Destructor
	virtual ~em() {};


protected:
	int m_iterations;							///< Number of iterations
	double m_loglikelihood;						///< Log likelihood
	double m_epsilon;							///< Small constant (1e-6)
	bool m_debug;								///< Debug mode
	OrderMethod m_order_method;					///< Variable ordering method
	InferMethod m_infer_method;					///< Inference method (for the E-step)
	int m_init_method;							///< Factor initialization method
	graphical_model m_gmo;						///< Original graphical model
	cte m_infer;								///< Inference engine (CTE)
	std::vector<std::vector<observation> > m_dataset;	///< Training dataset (-1 is a missing value)

private:

	std::vector<std::vector<vindex> > m_families;	///< Families X,pa(X): m_families[x] = pa(x)
	std::vector<factor> m_counts;					///< Counts for each family X,pa(X)
	std::string m_properties;
};


} // end namespace


#endif /* IBM_MERLIN_CORE_EM_H_ */
