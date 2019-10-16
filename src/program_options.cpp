/*
 * program_options.cpp
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

/// \file program_options.cpp
/// \brief Program options definitions
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#include "program_options.h"

///
/// \brief Parse the command line arguments.
///
ProgramOptions* parseCommandLine(int argc, char** argv) {

	ProgramOptions* opt = new ProgramOptions;

	// executable name
	opt->executableName = argv[0];

	try {

		po::options_description desc("Valid options");
		desc.add_options()
			("input-file,f", po::value<std::string>(), "path to problem file (required)")
			("evidence-file,e", po::value<std::string>(), "path to evidence file (required)")
			("query-file,q", po::value<std::string>(), "path to query file file (optional)")
			("virtual-evidence-file,V", po::value<std::string>(), "path to virtual evidence file (optional)")
			("output-file,o", po::value<std::string>(), "path to output file (optional)")
			("dataset-file,d", po::value<std::string>(), "path to dataset file (optional)")
			("algorithm,a", po::value<std::string>(), "inference algorithm (required): bte, cte, wmb, ijgp, lbp, jglp, gibbs")
			("task,t", po::value<std::string>(), "inference task (use PR, MAR, MAP, MMAP)")
			("ibound,i", po::value<int>(), "mini-bucket i-bound")
			("time-limit,l", po::value<int>(), "time limit in seconds")
			("seed,s", po::value<int>(), "seed for the random number generator")
			("verbose,v", po::value<int>(), "specify verbosity level")
			("debug,d", "enable debug mode")
			("positive,p", "enable positive probability values (> 0)")
			("iterations,n", po::value<int>(), "number of iterations")
			("samples,m", po::value<int>(), "number of samples")
			("threshold,E", po::value<double>(), "threshold value (default 1e-6)")
			("alpha,A", po::value<double>(), "equivalent sample size (default 5.0)")
			("init-factors,F", po::value<std::string>(), "initialize the factors")
			("output-format,O", po::value<std::string>(), "output file format (required)")
			("help,h", "produces this help message");

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		// parse help
		if (vm.count("help")) {
			std::cout << std::endl << desc << std::endl;
			delete opt;
			exit(0);
		}

		// parse verbosity level
		if (vm.count("verbose")) {
			opt->verbose = vm["verbose"].as<int>();
		}

		// parse debug mode (switch)
		if (vm.count("debug")) {
			opt->debug = true;
		}

		// parse positive mode (switch)
		if (vm.count("positive")) {
			opt->positive = true;
		}

		// parse input file
		if (!vm.count("input-file")) {
			std::string err_msg("Input model file is required. ");
			err_msg += "Call with '" + std::string(argv[0]) + " --help' ";
			err_msg += "for full description of the command line arguments.";
			throw std::runtime_error(err_msg);
		} else {
			opt->modelFile = vm["input-file"].as<std::string>();
		}

		// parse the query variables file
		if (vm.count("query-file")) {
			opt->queryFile = vm["query-file"].as<std::string>();
		}

		// parse the evidence file
		if (vm.count("evidence-file")) {
			opt->evidenceFile = vm["evidence-file"].as<std::string>();
		}

		// parse the virtual evidence file
		if (vm.count("virtual-evidence-file")) {
			opt->virtualEvidenceFile = vm["virtual-evidence-file"].as<std::string>();
		}

		// parse output file
		if (vm.count("output-file")) {
			opt->outputFile = vm["output-file"].as<std::string>();
		}

		// parse the data file
		if (vm.count("dataset-file")) {
			opt->datasetFile = vm["dataset-file"].as<std::string>();
		}

		// parse inference task
		if (vm.count("task")) {
			std::string task = vm["task"].as<std::string>();
			if (task.compare("PR") == 0) {
				opt->task = MERLIN_TASK_PR;
			} else if (task.compare("MAR") == 0) {
				opt->task = MERLIN_TASK_MAR;
			} else if (task.compare("MAP") == 0) {
				opt->task = MERLIN_TASK_MAP;
			} else if (task.compare("MMAP") == 0) {
				opt->task = MERLIN_TASK_MMAP;
			} else if (task.compare("EM") == 0) {
				opt->task = MERLIN_TASK_EM;
			} else {
				std::string err_msg("Inference task ");
				err_msg += task + " is not supported.";
				throw std::runtime_error(err_msg);
			}
		}

		// parse algorithm
		if (vm.count("algorithm")) {
			std::string alg = vm["algorithm"].as<std::string>();
			if (alg.compare("bte") == 0) {
				opt->algorithm = MERLIN_ALGO_BTE;
			} else if (alg.compare("cte") == 0) {
				opt->algorithm = MERLIN_ALGO_CTE;
			} else if (alg.compare("ijgp") == 0) {
				opt->algorithm = MERLIN_ALGO_IJGP;
			} else if (alg.compare("jglp") == 0) {
				opt->algorithm = MERLIN_ALGO_JGLP;
			} else if (alg.compare("gibbs") == 0) {
				opt->algorithm = MERLIN_ALGO_GIBBS;
			} else if (alg.compare("lbp") == 0) {
				opt->algorithm = MERLIN_ALGO_LBP;
			} else if (alg.compare("aobb") == 0) {
				opt->algorithm = MERLIN_ALGO_AOBB;
			} else if (alg.compare("aobf") == 0) {
				opt->algorithm = MERLIN_ALGO_AOBF;
			} else if (alg.compare("rbfaoo") == 0) {
				opt->algorithm = MERLIN_ALGO_RBFAOO;
			} else if (alg.compare("wmb") == 0) {
				opt->algorithm = MERLIN_ALGO_WMB;
			} else {
				std::string err_msg("Algorithm ");
				err_msg += alg + " is not supported.";
				throw std::runtime_error(err_msg);
			}
		}

		// parse mini-bucket i-bound
		if (vm.count("ibound")) {
			opt->ibound = vm["ibound"].as<int>();
		}

		// parse the time limit
		if (vm.count("time-limit")) {
			opt->timeLimit = vm["time-limit"].as<int>();
		}

		// parse the random generator seed
		if (vm.count("seed")) {
			opt->seed = vm["seed"].as<int>();
		}

		// parse the number of iterations
		if (vm.count("iterations")) {
			opt->iterations = vm["iterations"].as<int>();
		}

		// parse the number of samples
		if (vm.count("samples")) {
			opt->samples = vm["samples"].as<int>();
		}

		// parse the epsilon threshold
		if (vm.count("threshold")) {
			opt->threshold = vm["threshold"].as<double>();
		}

		// parse the equivalent sample size
		if (vm.count("alpha")) {
			opt->alpha = vm["alpha"].as<double>();
		}

		// parse the factor initialization
		if (vm.count("init-factors")) {
			std::string str = vm["init-factors"].as<std::string>();
			if (str.compare("none") == 0) {
				opt->initFactors = MERLIN_INIT_NONE;
			} else if (str.compare("random") == 0) {
				opt->initFactors = MERLIN_INIT_RANDOM;
			} else if (str.compare("uniform") == 0) {
				opt->initFactors = MERLIN_INIT_UNIFORM;
			} else {
				std::string err_msg("Factor initialization method ");
				err_msg += str + " is not supported.";
				throw std::runtime_error(err_msg);
			}
		}

		// parse the output format
		if (vm.count("output-format")) {
			std::string format = vm["output-format"].as<std::string>();
			if (format.compare("uai") == 0) {
				opt->outputFormat = MERLIN_OUTPUT_UAI;
			} else if (format.compare("json") == 0) {
				opt->outputFormat = MERLIN_OUTPUT_JSON;
			} else {
				std::string err_msg("The output format ");
				err_msg += format + " is not supported";
				throw std::runtime_error(err_msg);
			}
		}
	} catch (std::exception& e) {
		std::cerr << e.what() << std::endl;
		delete opt;
		return NULL;
	}

	return opt;
}

