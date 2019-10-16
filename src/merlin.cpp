/*
 * merlin.cpp
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

// Merlin library core.
#include "jglp.h"
#include "ijgp.h"
#include "lbp.h"
#include "gibbs.h"
#include "wmb.h"
#include "bte.h"
#include "cte.h"
#include "em.h"
#include "util.h"

#include "merlin.h"


///
/// \brief Constructs the default Merlin engine.
///
Merlin::Merlin() {
	m_task = MERLIN_TASK_MAR;
	m_algorithm = MERLIN_ALGO_WMB;
	m_ibound = 4;
	m_iterations = 100;
	m_samples = 100;
	m_gmo = NULL;
	m_debug = false;
	m_useFiles = true; // default: use input files
	m_outputFormat = MERLIN_OUTPUT_UAI;
	m_ioTime = 0;
	m_positive = false;
	m_threshold = 1e-06;
}

///
/// \brief Destroys the Merlin engine.
///
Merlin::~Merlin() {
	clear();
}

///
/// \brief Clears the internal structures.
///
void Merlin::clear() {
	if (m_gmo != NULL) {
		//delete static_cast<merlin::graphical_model*>(m_gmo);
		delete m_gmo;
		m_gmo = NULL;
	}
}

///
/// \brief Set the inference algorithm.
/// \param alg 	The code associated with the algorithm.
///
void Merlin::set_algorithm(size_t alg) {
	m_algorithm = alg;
}

///
/// \brief Set the inference task.
/// \param task	The code associated with the task.
///
void Merlin::set_task(size_t task) {
	m_task = task;
}

///
/// \brief Set the i-bound.
/// \param ibound The value of the i-bound parameter.
///
void Merlin::set_ibound(size_t ibound) {
	m_ibound = ibound;
}

///
/// \brief Set the number of iterations.
/// \param iter	The number of iterations.
///
void Merlin::set_iterations(size_t iter) {
	m_iterations = iter;
}

///
/// \brief Set the number of samples.
/// \param s	The number of samples.
///
void Merlin::set_samples(size_t s) {
	m_samples = s;
}

///
/// \brief Set the debug flag.
///
void Merlin::set_debug(bool v) {
	m_debug = v;
}

///
/// \brief Set the positive flag.
///
void Merlin::set_positive(bool v) {
	m_positive = v;
}

///
/// \brief Set the threshold value.
///
void Merlin::set_threshold(double e) {
	m_threshold = e;
}

///
/// \brief Set the equivalent sample size value.
///
void Merlin::set_alpha(double a) {
	m_alpha = a;
}

///
/// \brief Set the factor initialization method.
///
void Merlin::set_init_factor_method(int m) {
	m_initFactors = m;
}

///
/// \brief Set the use files flag
///
void Merlin::set_use_files(bool f) {
	m_useFiles = f;
}

///
/// \brief Set the input file name.
/// \param f	The file name.
///
void Merlin::set_model_file(std::string f) {
	m_modelFile = f;
}

///
/// \brief Set the output file name.
/// \param f	The file name.
///
void Merlin::set_output_file(std::string f) {
	m_outputFile = f;
}

///
/// \brief Set the evidence file name.
/// \param f	The file name.
///
void Merlin::set_evidence_file(std::string f) {
	m_evidenceFile = f;
}

///
/// \brief Set the virtual evidence file name.
/// \param f	The file name.
///
void Merlin::set_virtual_evidence_file(std::string f) {
	m_virtualEvidenceFile = f;
}

///
/// \brief Set the query file name.
/// \param f	The file name.
///
void Merlin::set_query_file(std::string f) {
	m_queryFile = f;
}

///
/// \brief Set the dataset file name.
/// \param f	The file name.
///
void Merlin::set_dataset_file(std::string f) {
	m_datasetFile = f;
}

///
/// \brief Set the input model string.
/// \param s	The model.
///
void Merlin::set_model_string(std::string s) {
	m_modelString = s;
}

///
/// \brief Set the output string.
/// \param s	The output.
///
void Merlin::set_output_string(std::string s) {
	m_outputString = s;
}

///
/// \brief Set the evidence string.
/// \param s	The evidence.
///
void Merlin::set_evidence_string(std::string s) {
	m_evidenceString = s;
}

///
/// \brief Set the query string.
/// \param s	The query.
///
void Merlin::set_query_string(std::string s) {
	m_queryString = s;
}

///
/// \brief Set the dataset string.
/// \param s	The dataset.
///
void Merlin::set_dataset_string(std::string s) {
	m_datasetString = s;
}

///
/// \brief Set the output format.
///
void Merlin::set_output_format(int f) {
	m_outputFormat = f;
}

///
/// \brief Read the graphical model.
/// \param filename	The input file name.
///
bool Merlin::read_model(const char* filename) {
	try {

		// Read the graphical model
		m_filename = std::string(filename);
		std::ifstream is(filename);
		if (is.fail()) {
			std::string err_msg("Cannot open the input file: ");
			err_msg += std::string(filename);
			throw std::runtime_error(err_msg);
		}

		merlin::graphical_model gm;
		gm.read(is, m_positive); // throws a runtime_error in case of failure

		// Clear any previous graphical model
		clear();

		// Store the original graphical model (without evidence)
		m_gmo = gm.clone();

		return true;
	} catch (const std::runtime_error& e) {
		std::cerr << e.what() << std::endl;
		return false;
	}
}

///
/// \brief Read the graphical model.
/// \param filename	The input file name.
///
bool Merlin::read_model(std::string model) {
	try {

		// Read the graphical model
		int id = merlin::randi(12345678);
		std::stringstream ss;
		ss << "model-" << id << ".uai";
		m_filename = ss.str();
		std::istringstream is(model);
		if (is.fail()) {
			std::string err_msg("Cannot open the input model string");
			throw std::runtime_error(err_msg);
		}

		merlin::graphical_model gm;
		gm.read(is, m_positive); // throws a runtime_error in case of failure

		// Clear any previous graphical model
		clear();

		// Store the original graphical model (without evidence)
		m_gmo = gm.clone();

		return true;
	} catch (const std::runtime_error& e) {
		std::cerr << e.what() << std::endl;
		return false;
	}
}

///
/// \brief Read the evidence as variable-value pairs.
/// \param filename	The evidence file name.
///
bool Merlin::read_evidence(const char* filename) {
	try {

		// Open the evidence file
		std::ifstream in(filename);
		if (in.fail()) {
			std::string err_msg("Cannot open the evidence file: ");
			err_msg += filename;
			throw std::runtime_error(err_msg);
		}

		// Clear any previous evidence (m_evidence is a map)
		m_evidence.clear();

		// Read the evidence file
		int num_evid;
		in >> num_evid;
		for (int i = 0; i < num_evid; ++i) {
			vindex var;
			size_t val;
			in >> var >> val;
			m_evidence[var] = val;
		}

		// Close the evidence file
		in.close();

		return true;
	} catch (const std::runtime_error& e) {
		std::cerr << e.what() << std::endl;
		return false;
	}
}

///
/// \brief Read the virtual evidence as variable-likelihood pairs.
/// \param filename	The virtual evidence file name.
///
bool Merlin::read_virtual_evidence(const char* filename) {
	try {

		// Open the evidence file
		std::ifstream in(filename);
		if (in.fail()) {
			std::string err_msg("Cannot open the virtual evidence file: ");
			err_msg += filename;
			throw std::runtime_error(err_msg);
		}

		// Clear any previous virtual evidence
		m_virtualEvidence.clear();

		// Read the virtual evidence file
		int num_evid;
		in >> num_evid;
		for (int i = 0; i < num_evid; ++i) {
			vindex var;
			size_t dom;
			likelihood l;
			in >> var >> dom;
			l.resize(dom);
			for (size_t k = 0; k < dom; ++k) {
				in >> l[k];
			}

			m_virtualEvidence[var] = l;
		}

		// Close the evidence file
		in.close();

		return true;
	} catch (const std::runtime_error& e) {
		std::cerr << e.what() << std::endl;
		return false;
	}
}

///
/// \brief Read the evidence as variable-value pairs.
/// \param filename	The evidence string.
///
bool Merlin::read_evidence(std::string evidence) {
	try {

		// Open the evidence file
		std::istringstream in(evidence);
		if (in.fail()) {
			std::string err_msg("Cannot open the evidence string");
			throw std::runtime_error(err_msg);
		}

		// Clear any previous evidence
		m_evidence.clear();

		// Read the evidence file
		int num_evid;
		in >> num_evid;
		for (int i = 0; i < num_evid; ++i) {
			vindex var;
			size_t val;
			in >> var >> val;
			m_evidence[var] = val;
		}

		return true;
	} catch (const std::runtime_error& e) {
		std::cerr << e.what() << std::endl;
		return false;
	}
}

///
/// \brief Read the virtual evidence as variable-likelihood pairs.
/// \param filename	The virtual evidence file name.
///
bool Merlin::read_virtual_evidence(std::string evidence) {
	try {

		// Open the evidence file
		std::istringstream in(evidence);
		if (in.fail()) {
			std::string err_msg("Cannot open the virtual evidence string");
			throw std::runtime_error(err_msg);
		}

		// Clear any previous virtual evidence
		m_virtualEvidence.clear();

		// Read the virtual evidence file
		int num_evid;
		in >> num_evid;
		for (int i = 0; i < num_evid; ++i) {
			vindex var;
			size_t dom;
			likelihood l;
			in >> var >> dom;
			l.resize(dom);
			for (size_t k = 0; k < dom; ++k) {
				in >> l[k];
			}

			m_virtualEvidence[var] = l;
		}

		return true;
	} catch (const std::runtime_error& e) {
		std::cerr << e.what() << std::endl;
		return false;
	}
}

///
/// \brief Read the query variables (MMAP task only).
/// \param filename	The query file name.
///
bool Merlin::read_query(const char* filename) {
	try {

		// Open the query file
		std::ifstream in(filename);
		if (in.fail()) {
			throw std::runtime_error("Error while opening the query file.");
		}

		// Clear any previous query
		m_query.clear();

		// Read the query file
		int num_vars;
		in >> num_vars;
		for (int i = 0; i < num_vars; ++i) {
			vindex var;
			in >> var;
			m_query.push_back(var);
		}

		// Close the evidence file
		in.close();

		// Sort the query variables in ascending order
		std::sort(m_query.begin(), m_query.end());

		return true;
	} catch (const std::runtime_error& e) {
		std::cerr << e.what() << std::endl;
		return false;
	}
}

///
/// \brief Read the query variables (MMAP task only).
/// \param filename	The query string.
///
bool Merlin::read_query(std::string query) {
	try {

		// Open the query file
		std::istringstream in(query);
		if (in.fail()) {
			std::string err_msg("Cannot open the query string");
			throw std::runtime_error(err_msg);
		}

		// Clear any previous query
		m_query.clear();

		// Read the query file
		int num_vars;
		in >> num_vars;
		for (int i = 0; i < num_vars; ++i) {
			vindex var;
			in >> var;
			m_query.push_back(var);
		}

		return true;
	} catch (const std::runtime_error& e) {
		std::cerr << e.what() << std::endl;
		return false;
	}
}

///
/// \brief Read the training dataset as vectors of values (possibly missing).
/// \param filename	The dataset file name.
///
bool Merlin::read_dataset(const char* filename) {
	try {

		// Open the evidence file
		std::ifstream in(filename);
		if (in.fail()) {
			std::string err_msg("Cannot open the training dataset file: ");
			err_msg += filename;
			throw std::runtime_error(err_msg);
		}

		// Clear any previous training dataset
		m_dataset.clear();

		// Read the training dataset file line by line
		std::string line;
		while (std::getline(in, line)) {
			if (line.empty()) continue; // skip empty lines
			std::vector<std::string> tokens = merlin::split(line, ',');
			std::vector<merlin::observation> example;
			for (size_t i = 0; i < tokens.size(); ++i) {
				std::string token = tokens[i];
				if (token.compare("?") == 0) { // missing value
					example.push_back(merlin::observation(i));
				} else if (token.find("[") != token.npos) { // likelihood evidence
					size_t first = token.find('[');
					size_t last = token.find(']');
					std::string str = token.substr(first+1, last-first-1);
					std::vector<std::string> temp = merlin::split(str, ';');
					std::vector<double> likelihood;
					for (size_t k = 0; k < temp.size(); ++k) {
						double d = std::atof(temp[k].c_str());
						likelihood.push_back(d);
					}
					example.push_back(merlin::observation(i, likelihood));
				} else { // regular value
					int val = std::atoi(token.c_str());
					example.push_back(merlin::observation(i, val));
				}
			}

			m_dataset.push_back(example);
		}

		// Close the training dataset file
		in.close();

		return true;
	} catch (const std::runtime_error& e) {
		std::cerr << e.what() << std::endl;
		return false;
	}
}

///
/// \brief Read the training dataset as vectors of values (possibly missing).
///        Lines in the dataset string are separated by '|', and values by ','.
/// \param filename	The dataset string.
///
bool Merlin::read_dataset(std::string dataset) {
	try {

		// Open the evidence file
		std::istringstream in(dataset);
		if (in.fail()) {
			std::string err_msg("Cannot open the training dataset string");
			throw std::runtime_error(err_msg);
		}

		// Clear any previous training dataset
		m_dataset.clear();

		// Read the training dataset string line by line
		std::vector<std::string> lines = merlin::split(dataset, '|');
		for (size_t l = 0; l < lines.size(); ++l) {
			std::string line = lines[l];
			if (line.empty()) continue; // skip empty lines
			std::vector<std::string> tokens = merlin::split(line, ',');
			std::vector<merlin::observation> example;
			for (size_t i = 0; i < tokens.size(); ++i) {
				std::string token = tokens[i];
				if (token.compare("?") == 0) { // missing value
					example.push_back(merlin::observation(i));
				} else if (token.find("[") != token.npos) { // likelihood evidence
					size_t first = token.find('[');
					size_t last = token.find(']');
					std::string str = token.substr(first+1, last-first-1);
					std::vector<std::string> temp = merlin::split(str, ';');
					std::vector<double> likelihood;
					for (size_t k = 0; k < temp.size(); ++k) {
						double d = std::atof(temp[k].c_str());
						likelihood.push_back(d);
					}
					example.push_back(merlin::observation(i, likelihood));
				} else { // regular value
					int val = std::atoi(token.c_str());
					example.push_back(merlin::observation(i, val));
				}
			}

			m_dataset.push_back(example);
		}

		return true;
	} catch (const std::runtime_error& e) {
		std::cerr << e.what() << std::endl;
		return false;
	}
}



///
/// \brief Write the graphical model.
/// \param filename	The output file name.
///
bool Merlin::write_model(const char* filename) {
	try {

		// Open the output file stream
		std::ofstream os(filename);
		if (os.fail()) {
			std::string err_msg("Cannot open the output file: ");
			err_msg += filename;
			throw std::runtime_error(err_msg);
		}

		// Write the graphical model
		merlin::graphical_model gm;
		gm = *(static_cast<merlin::graphical_model*>(m_gmo));
		gm.write(os); // default is UAI format
		os.close();

		return true;
	} catch (const std::runtime_error& e) {
		std::cerr << e.what() << std::endl;
		return false;
	}
}

///
/// \brief Safety checks.
///
void Merlin::check() {
	if (m_task == MERLIN_TASK_PR) {
		if (m_algorithm != MERLIN_ALGO_WMB &&
			m_algorithm != MERLIN_ALGO_BTE) {
			std::string err_msg("For PR inference use WMB and BTE algorithms.");
			throw std::runtime_error(err_msg);
		}
	} else if (m_task == MERLIN_TASK_MAR) {
		if (m_algorithm != MERLIN_ALGO_WMB &&
				m_algorithm != MERLIN_ALGO_IJGP &&
				m_algorithm != MERLIN_ALGO_LBP &&
				m_algorithm != MERLIN_ALGO_GIBBS &&
				m_algorithm != MERLIN_ALGO_BTE &&
				m_algorithm != MERLIN_ALGO_CTE) {
			std::string err_msg("For MAR inference use WMB, IJGP, LBP, BTE, CTE and GIBBS algorithms.");
			throw std::runtime_error(err_msg);
		}
	} else if (m_task == MERLIN_TASK_MAP) {
		if (m_algorithm != MERLIN_ALGO_WMB &&
				m_algorithm != MERLIN_ALGO_JGLP &&
				m_algorithm != MERLIN_ALGO_IJGP &&
				m_algorithm != MERLIN_ALGO_GIBBS &&
				m_algorithm != MERLIN_ALGO_BTE) {
			std::string err_msg("For MAP inference use WMB, JGLP, IJGP and GIBBS algorithms.");
			throw std::runtime_error(err_msg);
		}
	} else if (m_task == MERLIN_TASK_MMAP) {
		if (m_algorithm != MERLIN_ALGO_WMB &&
			m_algorithm != MERLIN_ALGO_BTE) {
			std::string err_msg("For MMAP inference use WMB and BTE algorithms.");
			throw std::runtime_error(err_msg);
		}
	} else if (m_task == MERLIN_TASK_EM) {
		if (m_algorithm != MERLIN_ALGO_CTE &&
			m_algorithm != MERLIN_ALGO_BTE &&
			m_algorithm != MERLIN_ALGO_WMB) {
			std::string err_msg("For EM learning use WMB, CTE and BTE inference algorithms.");
			throw std::runtime_error(err_msg);
		}
	} else {
		std::string err_msg("Supported inference tasks are PR, MAR, MAP and MMAP.");
		throw std::runtime_error(err_msg);
	}
}

///
/// \brief Initialize the solver.
///
bool Merlin::init() {

	double timestamp = merlin::timeSystem();
	if (m_useFiles) { // input files
		if (!read_model(m_modelFile.c_str())) {
			return false;
		}
		if (m_evidenceFile.empty() == false) {
			if (!read_evidence(m_evidenceFile.c_str())) {
				return false;
			}
		}
		if (m_virtualEvidenceFile.empty() == false) {
			if (!read_virtual_evidence(m_virtualEvidenceFile.c_str())) {
				return false;
			}
		}
		if (m_queryFile.empty() == false) {
			if (!read_query(m_queryFile.c_str())) {
				return false;
			}
		}
		if (m_datasetFile.empty() == false) {
			if (!read_dataset(m_datasetFile.c_str())) {
				return false;
			}
		}
	} else { // input strings
		if (!read_model(m_modelString)) {
			return false;
		}
		if (m_evidenceString.empty()) {
			if (!read_evidence(m_evidenceString)) {
				return false;
			}
		}
		if (m_virtualEvidenceString.empty()) {
			if (!read_virtual_evidence(m_virtualEvidenceString)) {
				return false;
			}
		}
		if (!m_queryString.empty()) {
			if (!read_query(m_queryString)) {
				return false;
			}
		}
		if (!m_datasetString.empty()) {
			if (!read_dataset(m_datasetString)) {
				return false;
			}
		}
	}

	m_ioTime = (merlin::timeSystem() - timestamp);
	return true;
}


///
/// \brief Solve the inference task given current evidence.
///
int Merlin::run() {

	try {

		// Safety checks
		check();

		// Prologue
		std::cout << VERSIONINFO << std::endl << COPYRIGHT << std::endl;
		std::cout << "[MERLIN] Initialize Merlin engine ..." << std::endl;
		std::cout << "[MERLIN] + tasks supported  : PR, MAR, MAP, MMAP, EM" << std::endl;

		// Initialize the graphical model
		merlin::graphical_model gm;
		std::vector<merlin::factor> fs, nfs;
		std::map<vindex, vindex> old2new;
		std::set<vindex> dummies;
		gm = *(static_cast<merlin::graphical_model*>(m_gmo)); // original model

		// Assert the plain evidence only (i.e., remove observed nodes from the graph)
		bool plainEvidence = (m_evidence.empty() ? false : true);
		bool virtualEvidence = (m_virtualEvidence.empty() ? false : true);

		if ( plainEvidence && !virtualEvidence) {
			fs = gm.assert_evidence( m_evidence, old2new );
		} else if (!plainEvidence && !virtualEvidence) {
			fs = gm.get_factors();
			for (size_t v = 0; v < gm.nvar(); ++v) {
				old2new[v] = v;
			}
		} else {

			// Create new factors P(U|X) for each virtual evidence variable
			vindex idx = gm.nvar();
			std::map<vindex, likelihood>::iterator mi = m_virtualEvidence.begin();
			for (; mi != m_virtualEvidence.end(); ++mi) {
				vindex x = mi->first;		// virtual evidence variable
				likelihood l = mi->second;	// likelihood

				// Check if virtual evidence variable is also plain evidence
				if (m_evidence.find(x) != m_evidence.end()) {
					std::stringstream err_msg;
					err_msg << "Variable " << x << " cannot be both virtual and regular evidence.";
					throw std::runtime_error(err_msg.str());
				}

				merlin::variable xvar = gm.var(x);
				merlin::variable uvar(idx++, 2);
				merlin::variable_set vs;
				vs |= xvar;
				vs |= uvar;
				merlin::factor f(vs, 0.0);
				f.set_child(uvar.label());
				m_evidence[uvar.label()] = 0; // first value of the U variable
				dummies.insert(uvar.label()); // keep track of dummy variables

				for (size_t k = 0; k < l.size(); ++k) {
					f.set(k, l[k]);
					f.set(k + xvar.states(), 1.0 - l[k]);
				}

				nfs.push_back(f); // extra factors
			}

			// Extend the model with the new factors
			for (size_t i = 0; i < nfs.size(); ++i) {
				gm.add_factor(nfs[i]);
			}

			// Assert the evidence
			fs = gm.assert_evidence( m_evidence, old2new );
		}

		// Set the output file
		if (m_outputFile.empty()) {
			size_t found = m_filename.find_last_of("/");
			std::string prob_name = (found != std::string::npos) ?
					m_filename.substr(found + 1) : m_filename;
			m_outputFile = "./" + prob_name;
		}

		// Setup the solver to run
		if ( m_task == MERLIN_TASK_PR ) { // PR inference

			// Set the output format
			m_outputFile += ".PR";
			if (m_outputFormat == MERLIN_OUTPUT_JSON) {
				m_outputFile += ".json";
			}

			// Open the output stream
			std::ofstream out(m_outputFile.c_str());
			if (out.fail()) {
				std::string err_msg("Cannot open output file: ");
				err_msg += m_outputFile;
				throw std::runtime_error(err_msg);
			}

			if (m_algorithm == MERLIN_ALGO_WMB) {
				merlin::wmb s(fs);
				std::ostringstream oss;
				oss << "iBound=" << m_ibound << ","
					<< "Order=MinFill" << ","
					<< "OrderIter=100" << ","
					<< "Iter=" << m_iterations << ","
					<< "Task=PR,Debug=" << (m_debug ? 1 : 0);
				s.set_properties(oss.str());
				s.run();
				s.write_solution(out, m_evidence, old2new, gm, dummies, m_outputFormat);
			} else if (m_algorithm == MERLIN_ALGO_BTE) {
				merlin::bte s(fs);
				std::ostringstream oss;
				oss << "Order=MinFill" << ","
					<< "Task=PR";
				s.set_properties(oss.str());
				s.run();
				s.write_solution(out, m_evidence, old2new, gm, dummies, m_outputFormat);
			}

			out.close();
		} else if ( m_task == MERLIN_TASK_MAR ) { // MAR inference

			// Set the output format
			m_outputFile += ".MAR";
			if (m_outputFormat == MERLIN_OUTPUT_JSON) {
				m_outputFile += ".json";
			}

			// Open the output stream
			std::ofstream out(m_outputFile.c_str());
			if (out.fail()) {
				std::string err_msg("Cannot open output file: ");
				err_msg += m_outputFile;
				throw std::runtime_error(err_msg);
			}

			if (m_algorithm == MERLIN_ALGO_WMB) {
				merlin::wmb s(fs);
				std::ostringstream oss;
				oss << "iBound=" << m_ibound << ","
					<< "Order=MinFill" << ","
					<< "OrderIter=100" << ","
					<< "Iter=" << m_iterations << ","
					<< "Task=MAR,Debug=" << (m_debug ? 1 : 0);
				s.set_properties(oss.str());
				s.run();
				s.write_solution(out, m_evidence, old2new, gm, dummies, m_outputFormat);
			} else if (m_algorithm == MERLIN_ALGO_IJGP) {
				merlin::ijgp s(fs);
				std::ostringstream oss;
				oss << "iBound=" << m_ibound << ","
					<< "Order=MinFill" << ","
					<< "Iter=" << m_iterations << ","
					<< "Task=MAR,Debug=" << (m_debug ? 1 : 0);
				s.set_properties(oss.str());
				s.run();
				s.write_solution(out, m_evidence, old2new, gm, dummies, m_outputFormat);
			} else if (m_algorithm == MERLIN_ALGO_LBP) {
				merlin::lbp s(fs);
				std::ostringstream oss;
				oss << "Schedule=Fixed,Distance=HPM,StopIter="
					<< m_iterations << ",StopObj=-1,StopMsg=-1,Debug=0";
				s.set_properties(oss.str());
				s.run();
				s.write_solution(out, m_evidence, old2new, gm, dummies, m_outputFormat);
			} else if (m_algorithm == MERLIN_ALGO_GIBBS) {
				merlin::gibbs s(fs);
				std::ostringstream oss;
				oss << "Task=MAR,TempMin=1.0,TempMax=1.0" << ","
					<< "Iter=" << m_iterations << ","
					<< "Samples=" << m_samples << ","
					<< "Debug=" << (m_debug ? 1 : 0);
				s.set_properties(oss.str());
				s.run();
				s.write_solution(out, m_evidence, old2new, gm, dummies, m_outputFormat);
			} else if (m_algorithm == MERLIN_ALGO_BTE) {
				merlin::bte s(fs);
				std::ostringstream oss;
				oss << "Order=MinFill" << ","
					<< "Task=MAR,Debug=" << (m_debug ? 1 : 0);
				s.set_properties(oss.str());
				s.run();
				std::ofstream out(m_outputFile.c_str());
				s.write_solution(out, m_evidence, old2new, gm, dummies, m_outputFormat);
			} else if (m_algorithm == MERLIN_ALGO_CTE) {
				merlin::cte s(fs);
				std::ostringstream oss;
				oss << "Order=MinFill" << ","
					<< "Task=MAR,Debug=" << (m_debug ? 1 : 0);
				s.set_properties(oss.str());
				std::vector<size_t> qvars;
				for (size_t i = 0; i < m_query.size(); ++i) {
					vindex var = m_query[i];
					vindex nvar = old2new.at(var);
					qvars.push_back(nvar); // use the new index of the MAP vars
				}
				s.set_query(qvars);
				s.run();
				std::ofstream out(m_outputFile.c_str());
				s.write_solution(out, m_evidence, old2new, gm, dummies, m_outputFormat);
			}

			out.close();
		} else if ( m_task == MERLIN_TASK_MAP ) { // MAP inference

			// Set the output format
			m_outputFile += ".MAP";
			if (m_outputFormat == MERLIN_OUTPUT_JSON) {
				m_outputFile += ".json";
			}

			// Open the output stream
			std::ofstream out(m_outputFile.c_str());
			if (out.fail()) {
				std::string err_msg("Cannot open output file: ");
				err_msg += m_outputFile;
				throw std::runtime_error(err_msg);
			}

			if (m_algorithm == MERLIN_ALGO_WMB) {
				merlin::wmb s(fs);
				std::ostringstream oss;
				oss << "iBound=" << m_ibound << ","
					<< "Order=MinFill" << ","
					<< "OrderIter=100" << ","
					<< "Iter=" << m_iterations << ","
					<< "Task=MAP";
				s.set_properties(oss.str());
				std::vector<vindex> qvars;
				for (size_t i = 0; i < gm.nvar(); ++i) {
					if (m_evidence.find(i) == m_evidence.end()) {
						size_t nvar = old2new.at(i);
						qvars.push_back(nvar); // use the new index of the MAP vars
					}
				}
				s.set_query(qvars);
				s.run();
				s.write_solution(out, m_evidence, old2new, gm, dummies, m_outputFormat);
			} else if (m_algorithm == MERLIN_ALGO_JGLP) {
				merlin::jglp s(fs);
				std::ostringstream oss;
				oss << "iBound=" << m_ibound << ","
					<< "Order=MinFill" << ","
					<< "Iter=" << m_iterations;
				s.set_properties(oss.str());
				s.run();
				s.write_solution(out, m_evidence, old2new, gm, dummies, m_outputFormat);
			} else if (m_algorithm == MERLIN_ALGO_IJGP) {
				merlin::ijgp s(fs);
				std::ostringstream oss;
				oss << "iBound=" << m_ibound << ","
					<< "Order=MinFill" << ","
					<< "Iter=" << m_iterations << ","
					<< "Task=MAP";
				s.set_properties(oss.str());
				s.run();
				s.write_solution(out, m_evidence, old2new, gm, dummies, m_outputFormat);
			} else if (m_algorithm == MERLIN_ALGO_GIBBS) {
				merlin::gibbs s(fs);
				std::ostringstream oss;
				oss << "TempMin=1.0,TempMax=1.0,Best=0,Beliefs=1" << ","
					<< "nIter=" << m_iterations << ","
					<< "nSamples=" << m_samples;
				s.set_properties(oss.str());
				s.run();
				s.write_solution(out, m_evidence, old2new, gm, dummies, m_outputFormat);
			} else if (m_algorithm == MERLIN_ALGO_BTE) {
				merlin::bte s(fs);
				std::ostringstream oss;
				oss << "Order=MinFill" << ","
					<< "Task=MAP";
				s.set_properties(oss.str());
				std::vector<vindex> qvars;
				for (size_t i = 0; i < gm.nvar(); ++i) {
					if (m_evidence.find(i) == m_evidence.end()) {
						size_t nvar = old2new.at(i);
						qvars.push_back(nvar); // use the new index of the MAP vars
					}
				}
				s.set_query(qvars);
				s.run();
				s.write_solution(out, m_evidence, old2new, gm, dummies, m_outputFormat);
			}

			out.close();
			// follow with search-based AOBB, AOBF, RBFAOO
		} else if ( m_task == MERLIN_TASK_MMAP ) { // MMAP inference

			// Set the output format
			m_outputFile += ".MMAP";
			if (m_outputFormat == MERLIN_OUTPUT_JSON) {
				m_outputFile += ".json";
			}

			// Open the output stream
			std::ofstream out(m_outputFile.c_str());
			if (out.fail()) {
				std::string err_msg("Cannot open output file: ");
				err_msg += m_outputFile;
				throw std::runtime_error(err_msg);
			}

			if (m_algorithm == MERLIN_ALGO_WMB) {
				merlin::wmb s(fs);
				std::ostringstream oss;
				oss << "iBound=" << m_ibound << ","
					<< "Order=MinFill" << ","
					<< "OrderIter=100" << ","
					<< "Iter=" << m_iterations << ","
					<< "Task=MMAP";
				s.set_properties(oss.str());
				std::vector<size_t> qvars;
				for (size_t i = 0; i < m_query.size(); ++i) {
					vindex var = m_query[i];
					vindex nvar = old2new.at(var);
					qvars.push_back(nvar); // use the new index of the MAP vars
				}
				s.set_query(qvars);
				s.run();
				s.write_solution(out, m_evidence, old2new, gm, dummies, m_outputFormat);
			} else if (m_algorithm == MERLIN_ALGO_BTE) {
				merlin::bte s(fs);
				std::ostringstream oss;
				oss << "Order=MinFill" << ","
					<< "Task=MMAP";
				s.set_properties(oss.str());
				std::vector<size_t> qvars;
				for (size_t i = 0; i < m_query.size(); ++i) {
					vindex var = m_query[i];
					vindex nvar = old2new.at(var);
					qvars.push_back(nvar); // use the new index of the MAP vars
				}
				s.set_query(qvars);
				s.run();
				s.write_solution(out, m_evidence, old2new, gm, dummies, m_outputFormat);
			}

			out.close();
			// follow with search-based AOBB, AOBF, RBFAOO
		} else if ( m_task == MERLIN_TASK_EM ) { // EM parameter learning

			merlin::em s(gm);

			// Set the output
			m_outputFile += ".EM";
			std::ofstream out(m_outputFile.c_str());
			if (out.fail()) {
				std::string err_msg("Cannot open output file: ");
				err_msg += m_outputFile;
				throw std::runtime_error(err_msg);
			}

			std::string initMethod("None");
			if (m_initFactors == MERLIN_INIT_UNIFORM) {
				initMethod = "Uniform";
			} else if (m_initFactors == MERLIN_INIT_RANDOM) {
				initMethod = "Random";
			} else {
				initMethod = "None";
			}

			std::ostringstream oss;
			oss << "Order=MinFill" << ","
				<< "Infer=CTE" << ","
				<< "Iter=" << m_iterations << ","
				<< "Debug=" << (m_debug ? "1" : "0") << ","
				<< "Threshold=" << m_threshold << ","
				<< "Init=" << initMethod;
			s.set_properties(oss.str());
			s.set_dataset(m_dataset);
			s.run();
			s.write_solution(out, gm);

			out.close();
		}

		std::cout << "[MERLIN] I/O time is " << std::fixed
				<< std::setprecision(MERLIN_PRECISION) << m_ioTime
				<< " seconds" << std::endl;

		return EXIT_SUCCESS;

	} catch (const std::runtime_error& e) {
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
	}
}


