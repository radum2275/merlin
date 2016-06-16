/*
 * run.cpp
 *
 *  Created on: 15 May 2015
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


#include "merlin.h"

#include "factor.h"
#include "algorithm.h"
#include "wmb.h"
#include "set.h"

///
/// \brief Solve an inference task given some evidence.
///
int run(const char* inputFile,
		const char* evidenceFile,
		const char* queryFile,
		const char* outputFile,
		const char* task,
		const unsigned int ibound,
		const unsigned int iterations) {

	try {

		// Safety checks
		assert(inputFile != NULL);
		assert(evidenceFile != NULL);
		assert(queryFile != NULL);
		assert(outputFile != NULL);

		// Initialize the graphical model
		merlin::graphical_model gm;
		std::vector<merlin::factor> fs;
		std::map<size_t, size_t> evidence, old2new;
		gm.read( inputFile );
		if ( strlen(evidenceFile) > 0) {
			fs = gm.assert_evidence( evidenceFile, evidence, old2new );
		} else {
			fs = gm.get_factors();
			for (size_t v = 0; v < gm.nvar(); ++v) old2new[v] = v;
		}

		// Setup the solver to run
		if (strcmp(task, "PR") == 0) {
			merlin::wmb s(fs);
			std::ostringstream oss;
			oss << "iBound=" << ibound << ","
				<< "Order=MinFill" << ","
				<< "Iter=" << iterations << ","
				<< "Task=PR";
			s.set_properties(oss.str());
			s.run();
			s.write_solution(outputFile, evidence, old2new, gm);
		} else if (strcmp(task, "MAR") == 0) {
			merlin::wmb s(fs);
			std::ostringstream oss;
			oss << "iBound=" << ibound << ","
				<< "Order=MinFill" << ","
				<< "Iter=" << iterations << ","
				<< "Task=MAR";
			s.set_properties(oss.str());
			s.run();
			s.write_solution(outputFile, evidence, old2new, gm);

		} else if (strcmp(task, "MAP") == 0) {
			merlin::wmb s(fs);
			std::ostringstream oss;
			oss << "iBound=" << ibound << ","
				<< "Order=MinFill" << ","
				<< "Iter=" << iterations << ","
				<< "Task=MAP";
			s.set_properties(oss.str());
			std::vector<size_t> qvars;
			for (size_t i = 0; i < gm.nvar(); ++i) {
				if (evidence.find(i) == evidence.end()) {
					size_t nvar = old2new.at(i);
					qvars.push_back(nvar); // use the new index of the MAP vars
				}
			}
			s.set_query(qvars);
			s.run();
			s.write_solution(outputFile, evidence, old2new, gm);

		} else if (strcmp(task, "MMAP") == 0) {
			merlin::wmb s(fs);
			std::ostringstream oss;
			oss << "iBound=" << ibound << ","
				<< "Order=MinFill" << ","
				<< "Iter=" << iterations << ","
				<< "Task=MMAP";
			s.set_properties(oss.str());
			std::vector<size_t> qvars;
			std::ifstream in(queryFile);
			if (in.fail()) {
				throw std::runtime_error("Error while opening query file (MMAP)");
			}

			// read the query variables
			size_t nvars;
			in >> nvars;
			for (size_t i = 0; i < nvars; ++i) {
				size_t var;
				in >> var;
				size_t nvar = old2new.at(var);
				qvars.push_back(nvar); // use the new index of the MAP vars
			}
			s.set_query(qvars);
			s.run();
			s.write_solution(outputFile, evidence, old2new, gm);

		} else {
			throw std::runtime_error("Unknown inference task. Use PR, MAR, MAP, MMAP.");
		}

		return 0; // success
	} catch(std::exception& e) {
		std::cerr << e.what() << std::endl;
		return 1; // failure
	}
}
