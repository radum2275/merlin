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
#include <iostream>
#include "merlin.h"

// Debugging only
void demo_debug() {

	// Init parameters
	unsigned int ibound = 2;
	unsigned int iterations = 10;
	unsigned int samples = 1000;
	const char* inputFile = "/home/radu/git/merlin/test/mrf.wcnf.uai";
	const char* evidenceFile = "/home/radu/git/merlin/example/simple5.evid";
	const char* queryFile = "/home/radu/git/merlin/example/simple5.map";
	const char* outputFile = "/home/radu/git/merlin/example/simple5.out";

//	const char* inputFile = "/home/radu/Downloads/chain.xmlbif.UAI";
//	const char* evidenceFile = "/home/radu/Downloads/chain.xmlbif.UAI.EVID";
//	const char* queryFile = "/home/radu/Downloads/chain.xmlbif.UAI.QUERY";

	// MMAP task
//	run(inputFile, evidenceFile, queryFile, outputFileMMAP, "MAR", ibound, iterations);

	// Initialize the Merlin engine
	Merlin eng;
	eng.set_param_ibound(ibound);
	eng.set_param_iterations(iterations);
	eng.set_param_samples(samples);
	eng.read_model(inputFile);
	eng.read_evidence(evidenceFile);
//	eng.read_query(queryFile);
	eng.set_task(MERLIN_TASK_MAR);
	eng.set_algorithm(MERLIN_ALGO_GIBBS);
	eng.run();
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

// Demo: convert a weigted cnf file into a factor graph (UAI format)
void demo_wcnf2uai(const char* file_name) {
	Merlin eng;

	std::string f = std::string(file_name) + ".uai";

	eng.read_model(file_name, MERLIN_INPUT_WCNF);
	std::cout << "Read wcnf file: " << file_name << std::endl;
	eng.write_model(f.c_str());
	std::cout << "Wrote uai file: " << f << std::endl;
}


// Main
int main(int argc, char** argv) {

	// Call the 'run' function
	//demo_run();

	// Call Merlin API
	//demo_api();

	// Call the 'debug' function
	demo_debug();

//	demo_convert(argv[1]);

//	demo_fileformat();

//	demo_merlin(argv[1]);

//	demo_wcnf2uai(argv[1]);

//	demo_wcnf2net(argv[1]);

	return 0;
}


