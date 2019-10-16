/*
 * program_options.h
 *
 *  Created on: 10 September 2018
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

/// \file program_options.h
/// \brief Program options definitions
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_MERLIN_PROGRAM_OPTIONS_H_
#define IBM_MERLIN_PROGRAM_OPTIONS_H_

#include <boost/program_options.hpp>
#include <iostream>
#include <sstream>

#include "base.h"

namespace po = boost::program_options;

struct ProgramOptions {

	std::string executableName;		///< Program name
	double timeLimit;				///< Time limit (in seconds)
	double memoryLimit;				///< Memory limit (in Gigs)
	int ibound;						///< Mini-bucket i-bound
	int algorithm;					///< Algorithm
	int task;						///< Inference task
	std::string modelFile;			///< Model file
	std::string evidenceFile;		///< Evidence file
	std::string queryFile;			///< Query file (MAP variables)
	std::string outputFile;			///< Output file
	std::string datasetFile;		///< Training dataset for learning
	std::string virtualEvidenceFile;///< Virtual evidence file
	int seed; 						///< Random number generator seed
	bool debug;						///< Debug mode
	int verbose;					///< Verbosity level (1=low, 2=medium, 3=high)
	int iterations;					///< Iterations for WMB/IJGP/JGLP
	int samples;					///< Number of samples
	int outputFormat;				///< Output format
	bool positive;					///< Force positive probabilities (>0)
	double threshold;				///< Tolerance threshold value (default 1e-06)
	double alpha;					///< Equivalent sample size (for Bayesian parameter estimation)
	int initFactors;					///< Initialize the CPTs (for EM learning)

public:

	// default constructor
	ProgramOptions();
};

ProgramOptions* parseCommandLine(int argc, char** argv);

inline ProgramOptions::ProgramOptions() :
		timeLimit(MERLIN_UNKNOWN),
		memoryLimit(80),
		ibound(2),
		algorithm(MERLIN_UNKNOWN),
		task(MERLIN_UNKNOWN),
		seed(12345678),
		debug(false),
		verbose(0),
		iterations(10),
		samples(1000),
		outputFormat(MERLIN_OUTPUT_UAI),
		positive(false),
		threshold(1e-6),
		alpha(5.0),
		initFactors(MERLIN_INIT_UNIFORM) {};

#endif /* IBM_MERLIN_PROGRAM_OPTIONS_H_ */
