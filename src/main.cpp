/*
 * merlin.cpp
 *
 *  Created on: 20 August 2018
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

#include <stdlib.h>
#include <string.h>
#include <iostream>

#include "merlin.h"
#include "program_options.h"

std::string fileToString(std::string filename) {
	std::ostringstream oss(std::ios::out | std::ios::binary); // *** binary
	std::ifstream in(filename.c_str());

	std::string line;
	while (std::getline(in, line)) {
		oss << line << "\r\n";
	}

	return oss.str();
}

int main(int argc, char** argv) {

	ProgramOptions* opt = parseCommandLine(argc, argv);
	if (!opt) {
		std::cerr << "Invalid command line arguments." << std::endl;
		return EXIT_FAILURE;
	}

	// Set the output format if not set
	if (opt->outputFormat < 0) {
		opt->outputFormat = MERLIN_OUTPUT_UAI;
	}

	// Setup Merlin engine
	Merlin eng;
	eng.set_use_files(true);
	eng.set_output_format(opt->outputFormat);
	eng.set_model_file(opt->modelFile);
	eng.set_evidence_file(opt->evidenceFile);
	eng.set_virtual_evidence_file(opt->virtualEvidenceFile);
	eng.set_output_file(opt->outputFile);
	eng.set_query_file(opt->queryFile);
	eng.set_dataset_file(opt->datasetFile);
	eng.set_task(opt->task);
	eng.set_algorithm(opt->algorithm);
	eng.set_ibound(opt->ibound);
	eng.set_iterations(opt->iterations);
	eng.set_samples(opt->samples);
	eng.set_debug(opt->debug);
	eng.set_positive(opt->positive);
	eng.set_threshold(opt->threshold);
	eng.set_alpha(opt->alpha);
	eng.set_init_factor_method(opt->initFactors);

	// Run the inference
	eng.init();
	int status = eng.run();

	delete opt;
	return status;
}

