
// Merlin library core.
#include "graphical_model.h"
#include "algorithm.h"
#include "jglp.h"
#include "ijgp.h"
#include "lbp.h"
#include "gibbs.h"
#include "wmb.h"
#include "set.h"

#include "merlin.h"


///
/// \brief Constructs the default Merlin engine.
///
Merlin::Merlin() {
	m_task = MERLIN_TASK_MAR;
	m_algorithm = MERLIN_ALGO_WMB;
	m_param_ibound = 4;
	m_param_iterations = 100;
	m_param_samples = 100;
	m_gmo = NULL;
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
		delete static_cast<merlin::graphical_model*>(m_gmo);
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
void Merlin::set_param_ibound(size_t ibound) {
	m_param_ibound = ibound;
}

///
/// \brief Set the number of iterations.
/// \param iter	The number of iterations.
///
void Merlin::set_param_iterations(size_t iter) {
	m_param_iterations = iter;
}

///
/// \brief Set the number of samples.
/// \param iter	The number of samples.
///
void Merlin::set_param_samples(size_t s) {
	m_param_samples = s;
}

///
/// \brief Read the graphical model.
/// \param f	The input file name.
///
bool Merlin::read_model(const char* f) {
	try {

		// Read the graphical model
		m_filename = std::string(f);
		merlin::graphical_model gm;
		gm.read(f); // throws a runtime_error in case of failure

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
/// \brief Assert evidence.
/// \param f	The evidence file name.
///
bool Merlin::read_evidence(const char* f) {
	try {

		// Open the evidence file
		std::ifstream in(f);
		if (in.fail()) {
			throw std::runtime_error("Error while opening the evidence file.");
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

		// Close the evidence file
		in.close();

		return true;
	} catch (const std::runtime_error& e) {
		std::cerr << e.what() << std::endl;
		return false;
	}
}

///
/// \brief Read the query variables (MMAP task only).
/// \param f	The query file name.
///
bool Merlin::read_query(const char* f) {
	try {

		// Open the query file
		std::ifstream in(f);
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

		return true;
	} catch (const std::runtime_error& e) {
		std::cerr << e.what() << std::endl;
		return false;
	}
}

///
/// \brief Write the graphical model.
/// \param f	The output file name.
///
bool Merlin::write_model(const char* f) {
	try {

		// Write the graphical model
		merlin::graphical_model gm;
		gm = *(static_cast<merlin::graphical_model*>(m_gmo));

		gm.write(f);

		return true;
	} catch (const std::runtime_error& e) {
		std::cerr << e.what() << std::endl;
		return false;
	}
}

///
/// \brief Write the graphical model.
/// \param f	The output file name.
///
bool Merlin::write_model(const char* f, int format) {
	try {

		// Write the graphical model
		merlin::graphical_model gm;
		gm = *(static_cast<merlin::graphical_model*>(m_gmo));

		gm.write(f, format);

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
		if (m_algorithm != MERLIN_ALGO_WMB) {
			throw std::runtime_error("Unsupported PR inference algorithm. Use WMB.");
		}
	} else if (m_task == MERLIN_TASK_MAR) {
		if (m_algorithm != MERLIN_ALGO_WMB &&
				m_algorithm != MERLIN_ALGO_IJGP &&
				m_algorithm != MERLIN_ALGO_LBP &&
				m_algorithm != MERLIN_ALGO_GIBBS) {
			throw std::runtime_error("Unsupported MAR inference algorithm. Use WMB, IJGP, LBP, GIBBS.");
		}
	} else if (m_task == MERLIN_TASK_MAP) {
		if (m_algorithm != MERLIN_ALGO_WMB &&
				m_algorithm != MERLIN_ALGO_JGLP &&
				m_algorithm != MERLIN_ALGO_IJGP &&
				m_algorithm != MERLIN_ALGO_GIBBS) {
			throw std::runtime_error("Unsupported MAP inference algorithm. Use WMB, JGLP, IJGP, GIBBS.");
		}
	} else if (m_task == MERLIN_TASK_MMAP) {
		if (m_algorithm != MERLIN_ALGO_WMB) {
			throw std::runtime_error("Unsupported MMAP inference algorithm. Use WMB.");
		}
	} else {
		throw std::runtime_error("Unsupported inference task. Use PR, MAR, MAP, MMAP");
	}
}

///
/// \brief Solve the inference task given current evidence.
///
int Merlin::run() {

	try {

		// Safety checks
		check();

		// Initialize the graphical model
		merlin::graphical_model gm;
		std::vector<merlin::factor> fs;
		std::map<vindex, vindex> old2new;
		gm = *(static_cast<merlin::graphical_model*>(m_gmo));

		// Assert the evidence
		if ( m_evidence.empty() == false ) {
			fs = gm.assert_evidence( m_evidence, old2new );
		} else {
			fs = gm.get_factors();
			for (size_t v = 0; v < gm.nvar(); ++v) {
				old2new[v] = v;
			}
		}

		// Setup the solver to run
		if ( m_task == MERLIN_TASK_PR ) {
			std::string output_file = m_filename + ".merlin.PR";
			merlin::wmb s(fs);
			std::ostringstream oss;
			oss << "iBound=" << m_param_ibound << ","
				<< "Order=MinFill" << ","
				<< "Iter=" << m_param_iterations << ","
				<< "Task=PR";
			s.set_properties(oss.str());
			s.run();
			s.write(output_file.c_str(), m_evidence, old2new, gm);
		} else if ( m_task == MERLIN_TASK_MAR ) {
			std::string output_file = m_filename + ".merlin.MAR";
			if (m_algorithm == MERLIN_ALGO_WMB) {
				merlin::wmb s(fs);
				std::ostringstream oss;
				oss << "iBound=" << m_param_ibound << ","
					<< "Order=MinFill" << ","
					<< "Iter=" << m_param_iterations << ","
					<< "Task=MAR";
				s.set_properties(oss.str());
				s.run();
				s.write(output_file.c_str(), m_evidence, old2new, gm);
			} else if (m_algorithm == MERLIN_ALGO_IJGP) {
				merlin::ijgp s(fs);
				std::ostringstream oss;
				oss << "iBound=" << m_param_ibound << ","
					<< "Order=MinFill" << ","
					<< "Iter=" << m_param_iterations << ","
					<< "Task=MAR";
				s.set_properties(oss.str());
				s.run();
				s.write(output_file.c_str(), m_evidence, old2new, gm);
			} else if (m_algorithm == MERLIN_ALGO_LBP) {
				merlin::lbp s(fs);
				std::ostringstream oss;
				oss << "Schedule=Priority,Distance=HPM,stopIter="
					<< m_param_iterations << ",stopObj=-1,stopMsg=-1";
				s.set_properties(oss.str());
				s.run();
				s.write(output_file.c_str(), m_evidence, old2new, gm);
			} else if (m_algorithm == MERLIN_ALGO_GIBBS) {
				merlin::gibbs s(fs);
				std::ostringstream oss;
				oss << "TempMin=1.0,TempMax=1.0,Best=1,Beliefs=0" << ","
					<< "nIter=" << m_param_iterations << ","
					<< "nSamples=" << m_param_samples;
				s.set_properties(oss.str());
				s.run();
				s.write(output_file.c_str(), m_evidence, old2new, gm);
			}
		} else if ( m_task == MERLIN_TASK_MAP ) {
			std::string output_file = m_filename + ".merlin.MAP";
			if (m_algorithm == MERLIN_ALGO_WMB) {
				merlin::wmb s(fs);
				std::ostringstream oss;
				oss << "iBound=" << m_param_ibound << ","
					<< "Order=MinFill" << ","
					<< "Iter=" << m_param_iterations << ","
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
				s.write(output_file.c_str(), m_evidence, old2new, gm);
			} else if (m_algorithm == MERLIN_ALGO_JGLP) {
				merlin::jglp s(fs);
				std::ostringstream oss;
				oss << "iBound=" << m_param_ibound << ","
					<< "Order=MinFill" << ","
					<< "Iter=" << m_param_iterations;
				s.set_properties(oss.str());
				s.run();
				s.write(output_file.c_str(), m_evidence, old2new, gm);
			} else if (m_algorithm == MERLIN_ALGO_IJGP) {
				merlin::ijgp s(fs);
				std::ostringstream oss;
				oss << "iBound=" << m_param_ibound << ","
					<< "Order=MinFill" << ","
					<< "Iter=" << m_param_iterations << ","
					<< "Task=MAP";
				s.set_properties(oss.str());
				s.run();
				s.write(output_file.c_str(), m_evidence, old2new, gm);
			} else if (m_algorithm == MERLIN_ALGO_GIBBS) {
				merlin::gibbs s(fs);
				std::ostringstream oss;
				oss << "TempMin=1.0,TempMax=1.0,Best=0,Beliefs=1" << ","
					<< "nIter=" << m_param_iterations << ","
					<< "nSamples=" << m_param_samples;
				s.set_properties(oss.str());
				s.run();
				s.write(output_file.c_str(), m_evidence, old2new, gm);
			} // follow with search-based AOBB, AOBF, RBFAOO
		} else if ( m_task == MERLIN_TASK_MMAP ) {
			std::string output_file = m_filename + ".merlin.MMAP";
			if (m_algorithm == MERLIN_ALGO_WMB) {
				merlin::wmb s(fs);
				std::ostringstream oss;
				oss << "iBound=" << m_param_ibound << ","
					<< "Order=MinFill" << ","
					<< "Iter=" << m_param_iterations << ","
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
				s.write(output_file.c_str(), m_evidence, old2new, gm);
			} // follow with search-based AOBB, AOBF, RBFAOO
		}

		return EXIT_SUCCESS;
	} catch (const std::runtime_error& e) {
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
	}
}


