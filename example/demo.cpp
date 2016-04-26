/*
 ============================================================================
 Name        : demo.cpp
 Author      : Radu Marinescu
 Version     :
 Copyright   : Copyright (c) IBM Corp. 2015
 Description : Uses shared library to print greeting
 To run the resulting executable the LD_LIBRARY_PATH must be
 set to ${project_loc}/merlin/.libs
 Alternatively, libtool creates a wrapper shell script in the
 build directory of this program which can be used to run it.
 Here the script will be called exampleProgram.
 ============================================================================
 */
#include <string>
#include "merlin.h"

// Debugging only
void demo_debug() {

	// Init parameters
	unsigned int ibound = 2;
	unsigned int iterations = 100;
	const char* inputFile = "/home/radu/git/merlin/example/smoking3.uai";
	const char* evidenceFile = "/home/radu/git/merlin/example/smoking3.evid";
	const char* queryFile = "/home/radu/git/merlin/example/smoking3.map";
	const char* outputFileMAR = "/home/radu/git/merlin/example/smoking3.MAR.out";
	const char* outputFileMAP = "/home/radu/git/merlin/example/smoking3.MAP.out";
	const char* outputFileMMAP = "/home/radu/git/merlin/example/smoking3.MMAP.out";

	// MMAP task
	run(inputFile, evidenceFile, queryFile, outputFileMMAP, "MAR", ibound, iterations);

}

// Demo the black-box run
void demo_run() {

	// Init parameters
	unsigned int ibound = 4;
	unsigned int iterations = 300;
	const char* inputFile = "pedigree1.uai";
	const char* evidenceFile = "pedigree1.evid";
	const char* queryFile = "pedigree1.map";
	const char* outputFileMAR = "pedigree1.MAR.out";
	const char* outputFileMAP = "pedigree1.MAP.out";
	const char* outputFileMMAP = "pedigree1.MMAP.out";

	// MAR task
	run(inputFile, evidenceFile, "", outputFileMAR, "MAR", ibound, iterations);

	// MAP task
	run(inputFile, evidenceFile, "", outputFileMAP, "MAP", ibound, iterations);

	// MMAP task
	run(inputFile, evidenceFile, queryFile, outputFileMMAP, "MMAP", ibound, iterations);

}

// Demo the Merlin API
void demo_api() {

	// Init parameters
	unsigned int ibound = 4;
	unsigned int iterations = 300;
	const char* model_file = "pedigree1.uai";
	const char* evid_file = "pedigree1.evid";
	const char* query_file = "pedigree1.map";


	// Initialize the Merlin engine
	Merlin eng;
	eng.set_param_ibound(8);
	eng.set_param_iterations(300);
	eng.read_model(model_file);
	eng.read_evidence(evid_file);

	// Solve a MAR task
	eng.set_task(MERLIN_TASK_MAR);
	eng.set_algorithm(MERLIN_ALGO_WMB);
	eng.run();

	// Solve a MAP task
	eng.set_task(MERLIN_TASK_MAP);
	eng.set_algorithm(MERLIN_ALGO_WMB);
	eng.run();

	// Solve a MMAP task
	eng.read_query(query_file);
	eng.set_task(MERLIN_TASK_MMAP);
	eng.run();
}


// Demo the Merlin API for Alchemy2 input
void demo_convert(const char* file_name) {

	// Init parameters
//	const char* model_file = "/home/radu/research/ibm/projects/csi/cpa/praline/results/Q1/A/mrf-alc2-le-logs.uai";
//	const char* model_file = "/home/radu/git/merlin/test/mrf.uai";
//	const char* evid_file = "/home/radu/research/ibm/projects/csi/cpa/praline/results/Q1/B/mrf.uai.evid";
//	const char* query_file = "/home/radu/research/ibm/projects/csi/cpa/praline/results/Q1/B/mrf.uai.map";
//	const char* outputFileMAR = "/home/radu/git/merlin/example/smoking3.MAR.out";
//	const char* outputFileMAP = "/home/radu/git/merlin/example/smoking3.MAP.out";
//	const char* outputFileMMAP = "/home/radu/git/merlin/example/smoking3.MMAP.out";
//	const char* out_file = "/home/radu/git/merlin/test/mrf-merlin.uai";
//	const char* out_file0 = "/home/radu/research/ibm/projects/csi/cpa/praline/results/Q1/A/mrf-nl.uai";
//	const char* out_file1 = "/home/radu/research/ibm/projects/csi/cpa/praline/results/Q1/A/mrf-en.uai";
//	const char* out_file2 = "/home/radu/research/ibm/projects/csi/cpa/praline/results/Q1/A/mrf-eu.uai";

	// Initialize the Merlin engine
	Merlin eng;
//	eng.set_param_ibound(2);
//	eng.set_param_iterations(100);
	eng.read_model(file_name);
	//eng.read_evidence(evid_file);
	std::string f0 = std::string(file_name) + ".uai";
	std::string f2 = std::string(file_name) + "-eu.uai";
	eng.write_model(f0.c_str(), 0);
	eng.write_model(f2.c_str(), 2);

	// Solve a MAR task
//	eng.set_task(MERLIN_TASK_MAR);
//	eng.set_algorithm(MERLIN_ALGO_WMB);
//	eng.run();

	// Solve a MAP task
//	eng.set_task(MERLIN_TASK_MAP);
//	eng.set_algorithm(MERLIN_ALGO_WMB);
//	eng.run();

	// Solve a MMAP task
//	eng.read_query(query_file);
//	eng.set_task(MERLIN_TASK_MMAP);
//	eng.run();
}

// Demo the Merlin API for Alchemy2 input
void demo_fileformat() {

	// Init parameters
//	const char* model_file = "/home/radu/research/ibm/projects/csi/cpa/praline/results/Q1/A/mrf-alc2-le-logs.uai";
//	const char* model_file = "/home/radu/git/merlin/test/mrf.uai";
//	const char* out_file = "/home/radu/git/merlin/test/mrf-merlin.uai";

	const char* model_file = "/home/radu/git/merlin/test/mrf-merlin.uai";
	const char* out_file = "/home/radu/git/merlin/test/mrf-merlin2.uai";

	// Initialize the Merlin engine
	Merlin eng;
	eng.read_model(model_file);
	eng.write_model(out_file);

}

// Demo the Merlin API for Alchemy2 input
void demo_merlin(const char* file_name) {

	// Initialize the Merlin engine
	Merlin eng;
	eng.set_param_ibound(12);
	eng.set_param_iterations(100);
	eng.read_model(file_name);

	// Solve a MAR task
	eng.set_task(MERLIN_TASK_MAR);
	eng.set_algorithm(MERLIN_ALGO_WMB);
	eng.run();

	// Solve a MAP task
//	eng.set_task(MERLIN_TASK_MAP);
//	eng.set_algorithm(MERLIN_ALGO_WMB);
//	eng.run();

	// Solve a MMAP task
//	eng.read_query(query_file);
//	eng.set_task(MERLIN_TASK_MMAP);
//	eng.run();
}


// Main
int main(int argc, char** argv) {

	// Call the 'run' function
	//demo_run();

	// Call Merlin API
	//demo_api();

	// Call the 'debug' function
	//demo_debug();

	demo_convert(argv[1]);

//	demo_fileformat();

//	demo_merlin(argv[1]);

	return 0;
}


