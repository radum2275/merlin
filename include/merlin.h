/*
 * merlin.h
 *
 *  Created on: 20 May 2015
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

/// \file merlin.h
/// \brief Merlin inference engine interface
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef __IBM_MERLIN_H_
#define __IBM_MERLIN_H_

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <string>
#include <memory>

#include "graphical_model.h"
#include "observation.h"

///
/// Merlin probabilistic inference engine.
///
class Merlin {
	typedef size_t vindex;
	typedef std::vector<double> likelihood;

protected:
	// Members:

	size_t m_task;						///< Inference task (PR, MAR, MAP, MMAP).
	size_t m_algorithm;					///< Inference algorithm
	size_t m_ibound;					///< Parameter: i-bound
	size_t m_iterations;				///< Parameter: number of iterations
	size_t m_samples;					///< Parameter: number of samples
	std::string m_modelFile;			///< Model file name
	std::string m_evidenceFile;			///< Evidence file name
	std::string m_virtualEvidenceFile;	///< Virtual evidence file name
	std::string m_queryFile;			///< Query file name
	std::string m_outputFile;			///< Output file name
	std::string m_datasetFile;			///< Training dataset file name
	std::string m_modelString;			///< Model string
	std::string m_evidenceString;		///< Evidence string
	std::string m_virtualEvidenceString;///< Virtual evidence string
	std::string m_queryString;			///< Query string
	std::string m_outputString;			///< Output string
	std::string m_datasetString;		///< Training dataset string
	bool m_debug;						///< Debug mode
	bool m_useFiles;					///< Use files or strings for IO
	int m_outputFormat;					///< Output format (UAI, JSON)
	bool m_positive;					///< Positive mode (probabilities > 0)
	double m_threshold;					///< Threshold value
	double m_alpha;						///< Equivalent sample size
	int m_initFactors;					///< Factor initialization method

private:
	// Local members:

	merlin::graphical_model* m_gmo;					///< Original graphical model.
	std::map<vindex, size_t> m_evidence;			///< Evidence as variable value pairs.
	std::map<vindex, likelihood> m_virtualEvidence;	///< Virtual evidence as variable likelihood pairs.
	std::vector<vindex> m_query;					///< Query variables for MMAP tasks.
	std::string m_filename;							///< Input model file name.
	std::vector<std::vector<merlin::observation> > m_dataset;		///< Training dataset.
	double m_ioTime;

	///
	/// \brief Clear the existing graphical model.
	///
	void clear();

	///
	/// \brief Perform safety checks.
	////
	void check();

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
	void set_ibound(size_t ibound);

	///
	/// \brief Set the number of iterations.
	/// \param iter	The number of iterations.
	///
	void set_iterations(size_t iter);

	///
	/// \brief Set the number of samples.
	/// \param s	The number of samples.
	///
	void set_samples(size_t s);

	///
	/// \brief Set the debug mode.
	/// \param d	The flag.
	///
	void set_debug(bool v);

	///
	/// \brief Set the positive mode.
	/// \param v	The flag.
	///
	void set_positive(bool v);

	///
	/// \brief Set the threshold value.
	/// \param e	The threshold.
	///
	void set_threshold(double e);

	///
	/// \brief Set the equivalent sample size.
	/// \param a	The equivalent sample size.
	///
	void set_alpha(double a);

	///
	/// \brief Set the factor initialization method.
	/// \param m	The initialization method.
	///
	void set_init_factor_method(int m);

	///
	/// \brief Set the input file name.
	/// \param f	The file name.
	///
	void set_model_file(std::string f);

	///
	/// \brief Set the output file name.
	/// \param f	The file name.
	///
	void set_output_file(std::string f);

	///
	/// \brief Set the evidence file name.
	/// \param f	The file name.
	///
	void set_evidence_file(std::string f);

	///
	/// \brief Set the virtual evidence file name.
	/// \param f	The file name.
	///
	void set_virtual_evidence_file(std::string f);

	///
	/// \brief Set the query file name.
	/// \param f	The file name.
	///
	void set_query_file(std::string f);

	///
	/// \brief Set the dataset file name.
	/// \param f	The file name.
	///
	void set_dataset_file(std::string f);

	///
	/// \brief Set the input model string.
	/// \param f	The model.
	///
	void set_model_string(std::string s);

	///
	/// \brief Set the output string.
	/// \param f	The output.
	///
	void set_output_string(std::string s);

	///
	/// \brief Set the evidence string.
	/// \param f	The evidence.
	///
	void set_evidence_string(std::string s);

	///
	/// \brief Set the query string.
	/// \param f	The query.
	///
	void set_query_string(std::string s);

	///
	/// \brief Set the dataset string.
	/// \param f	The dataset.
	///
	void set_dataset_string(std::string s);

	///
	/// \brief Set the flag indicating input files or strings.
	/// \param f	The flag.
	///
	void set_use_files(bool f);

	///
	/// \brief Set output format.
	/// \param f	The format.
	///
	void set_output_format(int f);

	///
	/// \brief Initialize the solver.
	///	\return *true* if succesful and *false* otherwise.
	///
	bool init();

	///
	/// \brief Solve the inference task given current evidence.
	///	\return 0 if succesful and 1 otherwise.
	///
	int run();


protected:

	///
	/// \brief Read the graphical model from a file in the UAI format.
	/// \param file_name	The input file name.
	///	\return *true* if successful and *false* otherwise.
	///
	bool read_model(const char* filename);

	///
	/// \brief Read the evidence.
	/// \param file_name	The evidence file name.
	///	\return *true* if successful and *false* otherwise.
	//
	bool read_evidence(const char* filename);

	///
	/// \brief Read the virtual evidence.
	/// \param file_name	The virtual evidence file name.
	///	\return *true* if successful and *false* otherwise.
	//
	bool read_virtual_evidence(const char* filename);

	///
	/// \brief Read the query variables (MMAP task and joint marginals).
	/// \param file_name	The query file name.
	///	\return *true* if successful and *false* otherwise.
	///
	bool read_query(const char* filename);

	///
	/// \brief Read the training dataset.
	/// \param file_name	The dataset file name.
	///	\return *true* if successful and *false* otherwise.
	///
	bool read_dataset(const char* filename);

	///
	/// \brief Read the graphical model from a file in the UAI format.
	/// \param model	The model string.
	///	\return *true* if successful and *false* otherwise.
	///
	bool read_model(std::string model);

	///
	/// \brief Read the evidence as variable-value pairs.
	/// \param evidence	The evidence string.
	///	\return *true* if successful and *false* otherwise.
	//
	bool read_evidence(std::string evidence);

	///
	/// \brief Read the virtual evidence as variable-likelihood pairs.
	/// \param evidence	The virtual evidence string.
	///	\return *true* if successful and *false* otherwise.
	//
	bool read_virtual_evidence(std::string evidence);

	///
	/// \brief Read the query variables (MMAP task only).
	/// \param query	The query string.
	///	\return *true* if successful and *false* otherwise.
	///
	bool read_query(std::string query);

	///
	/// \brief Read the training dataset.
	/// \param query	The dataset string.
	///	\return *true* if successful and *false* otherwise.
	///
	bool read_dataset(std::string dataset);

	///
	/// \brief Write the graphical model to a file in the specified format
	/// \param f		The output file name.
	/// \param format	The file format supported.
	///	\return *true* if successful and *false* otherwise.
	///
	bool write_model(const char* filename);

};

#endif /* __IBM_MERLIN_H_ */
