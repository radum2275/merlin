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

#include "merlin.h"

// Debugging only
void demo_debug() {

	// Init parameters
	unsigned int ibound = 2;
	unsigned int iterations = 300;
	const char* inputFile = "/home/radu/Downloads/chain.xmlbif.UAI.1";
	const char* evidenceFile = "/home/radu/Downloads/chain.xmlbif.UAI.EVID";
	const char* queryFile = "/home/radu/Downloads/chain.xmlbif.UAI.QUERY";
	const char* outputFileMAR = "pedigree1.MAR.out";
	const char* outputFileMAP = "pedigree1.MAP.out";
	const char* outputFileMMAP = "/home/radu/Downloads/chain.xmlbif.UAI.MMAP.out";

	// MMAP task
	run(inputFile, evidenceFile, queryFile, outputFileMMAP, "MMAP", ibound, iterations);

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
	eng.set_param_ibound(4);
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

// Main
int main(void) {

	// Call the 'run' function
	//demo_run();

	// Call Merlin API
	//demo_api();

	// Call the 'debug' function
	demo_debug();

	return 0;
}
