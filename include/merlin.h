/*
 * merlin.h
 *
 *  Created on: 20 May 2015
 *      Author: radu
 */

#ifndef __IBM_MERLIN_H_
#define __IBM_MERLIN_H_

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <string>
#include <memory>

///
/// Probabilistic inference algorithms.
///
#define MERLIN_ALGO_GIBBS 	1000		///< Gibbs sampling
#define MERLIN_ALGO_LBP		1001		///< Loopy BP
#define MERLIN_ALGO_IJGP	1002		///< Iterative join graph propagation
#define MERLIN_ALGO_JGLP	1003		///< Join graph linear programming
#define MERLIN_ALGO_WMB		1004		///< Weighted mini-buckets
#define MERLIN_ALGO_AOBB	1005		///< AND/OR branch and bound
#define MERLIN_ALGO_AOBF	1006		///< AO*
#define MERLIN_ALGO_RBFAOO	1007		///< Recursive best-first AND/OR search

///
/// Probabilistic inference tasks.
///
#define MERLIN_TASK_PR		10			///< Partition function (probability of evidence)
#define MERLIN_TASK_MAR		20			///< Posterior marginals (given evidence)
#define MERLIN_TASK_MAP		30			///< Maximum aposteriori (given evidence)
#define MERLIN_TASK_MMAP	40			///< Marginal MAP (given evidence)

///
/// Merlin probabilistic inference engine.
///
class Merlin {
	typedef size_t vindex;

protected:
	// Members:

	size_t m_task;					///< Inference task (PR, MAR, MAP, MMAP).
	size_t m_algorithm;				///< Inference algorithm
	size_t m_param_ibound;			///< Parameter: i-bound
	size_t m_param_iterations;		///< Parameter: iterations
	size_t m_param_samples;			///< Number of samples (Gibbs only)

private:
	// Local members:

	void* m_gmo;						///< Original graphical model.
	std::map<vindex, size_t> m_evidence;///< Evidence as variable value pairs.
	std::vector<vindex> m_query;		///< Query variables for MMAP tasks.
	std::string m_filename;				///< Input file name.
	void clear();						///< Clear the model.
	void check();						///< Perform safety checks

public:

	///
	/// \brief Constructs the default Merlin engine.
	///
	Merlin();

	///
	/// \brief Destroys the Merlin engine.
	///
	~Merlin();

	///
	/// \brief Set the inference algorithm.
	/// \param alg 	The code associated with the algorithm.
	///
	void set_algorithm(size_t alg);

	///
	/// \brief Set the inference task.
	/// \param task	The code associated with the task.
	///
	void set_task(size_t task);

	///
	/// \brief Set the i-bound.
	/// \param ibound The value of the i-bound parameter.
	///
	void set_param_ibound(size_t ibound);

	///
	/// \brief Set the number of iterations.
	/// \param iter	The number of iterations.
	///
	void set_param_iterations(size_t iter);

	///
	/// \brief Set the number of samples.
	/// \param s	The number of samples.
	///
	void set_param_samples(size_t s);

	///
	/// \brief Read the graphical model.
	/// \param f	The input file name.
	///	\return *true* if successful and *false* otherwise.
	///
	bool read_model(const char* f);

	///
	/// \brief Assert evidence.
	/// \param f	The evidence file name.
	///	\return *true* if successful and *false* otherwise.
	//
	bool read_evidence(const char* f);

	///
	/// \brief Read the query variables (MMAP task only).
	/// \param f	The query file name.
	///	\return *true* if successful and *false* otherwise.
	///
	bool read_query(const char* f);

	///
	/// \brief Solve the inference task given current evidence.
	///	\return 0 if succesful and 1 otherwise.
	///
	int run();

	///
	/// \brief Write the graphical model.
	/// \param f	The output file name.
	///	\return *true* if successful and *false* otherwise.
	///
	bool write_model(const char* f);
	bool write_model(const char* f, int format);

};


// Solves the inference task (black box)
extern "C" {
int run(const char* inputFile,
		const char* evidenceFile,
		const char* queryFile,
		const char* outputFile,
		const char* task,
		const unsigned int ibound = 2,
		const unsigned int = 100);
}

#endif /* __IBM_MERLIN_H_ */
