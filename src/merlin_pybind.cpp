#include <pybind11/pybind11.h>
#include <pybind11/attr.h>
#include <pybind11/stl.h>

#include "merlin.h"
#include "base.h"

namespace py = pybind11;

PYBIND11_MODULE(merlin, m) {
    m.doc() = "Merlin probabilistic inference engine binding. "
              "Provides an interface to the Merlin C++ library for various "
              "inference tasks on graphical models.";

    // Bind the Merlin class
    py::class_<Merlin>(m, "Merlin", "Merlin probabilistic inference engine.")
        .def(py::init<>(), "Constructs the default Merlin engine.")
        .def_property("task", &Merlin::get_task, &Merlin::set_task,
                      "Inference task (PR, MAR, MAP, MMAP).")
        .def_property("algorithm", &Merlin::get_algorithm, &Merlin::set_algorithm,
                      "Inference algorithm.")
        .def_property("ibound", &Merlin::get_ibound, &Merlin::set_ibound,
                      "Parameter: i-bound.")
        .def_property("iterations", &Merlin::get_iterations, &Merlin::set_iterations,
                      "Parameter: number of iterations.")
        .def_property("samples", &Merlin::get_samples, &Merlin::set_samples,
                      "Parameter: number of samples.")
        .def_property("model_file", &Merlin::get_model_file, &Merlin::set_model_file,
                      "Model file name.")
        .def_property("evidence_file", &Merlin::get_evidence_file, &Merlin::set_evidence_file,
                      "Evidence file name.")
        .def_property("virtual_evidence_file", &Merlin::get_virtual_evidence_file, &Merlin::set_virtual_evidence_file,
                      "Virtual evidence file name.")
        .def_property("query_file", &Merlin::get_query_file, &Merlin::set_query_file,
                      "Query file name.")
        .def_property("output_file", &Merlin::get_output_file, &Merlin::set_output_file,
                      "Output file name.")
        .def_property("dataset_file", &Merlin::get_dataset_file, &Merlin::set_dataset_file,
                      "Training dataset file name.")
        .def_property("model_string", &Merlin::get_model_string, &Merlin::set_model_string,
                      "Model string.")
        .def_property("evidence_string", &Merlin::get_evidence_string, &Merlin::set_evidence_string,
                      "Evidence string.")
        .def_property("virtual_evidence_string", &Merlin::get_virtual_evidence_string, &Merlin::set_virtual_evidence_string,
                      "Virtual evidence string.")
        .def_property("query_string", &Merlin::get_query_string, &Merlin::set_query_string,
                      "Query string.")
        .def_property("output_string", &Merlin::get_output_string, &Merlin::set_output_string,
                      "Output string.")
        .def_property("dataset_string", &Merlin::get_dataset_string, &Merlin::set_dataset_string,
                      "Training dataset string.")
        .def_property("debug", &Merlin::get_debug, &Merlin::set_debug,
                      "Debug mode flag.")
        .def_property("use_files", &Merlin::get_use_files, &Merlin::set_use_files,
                      "Flag: Use files or strings for IO.")
        .def_property("output_format", &Merlin::get_output_format, &Merlin::set_output_format,
                      "Output format (UAI, JSON).")
        .def_property("positive", &Merlin::get_positive, &Merlin::set_positive,
                      "Positive mode (probabilities > 0).")
        .def_property("threshold", &Merlin::get_threshold, &Merlin::set_threshold,
                      "Threshold value.")
        .def_property("alpha", &Merlin::get_alpha, &Merlin::set_alpha,
                      "Equivalent sample size.")
        .def_property("init_factors", &Merlin::get_init_factor_method, &Merlin::set_init_factor_method,
                      "Factor initialization method.")
        .def("init", &Merlin::init,
             "Initialize the solver.\n\n"
             ":return: True if successful, False otherwise.")
        .def("run", &Merlin::run,
             "Solve the inference task given current evidence.\n\n"
             ":return: 0 if successful, 1 otherwise.");

    //  --------- Enumerations ----------
    py::object enum_module = py::module_::import("enum");
    py::object PyIntEnum = enum_module.attr("IntEnum");

    py::dict algorithm_members;
    algorithm_members[py::str("GIBBS")] = py::cast(MERLIN_ALGO_GIBBS);
    algorithm_members[py::str("LBP")] = py::cast(MERLIN_ALGO_LBP);
    algorithm_members[py::str("IJGP")] = py::cast(MERLIN_ALGO_IJGP);
    algorithm_members[py::str("JGLP")] = py::cast(MERLIN_ALGO_JGLP);
    algorithm_members[py::str("WMB")] = py::cast(MERLIN_ALGO_WMB);
    algorithm_members[py::str("AOBB")] = py::cast(MERLIN_ALGO_AOBB);
    algorithm_members[py::str("AOBF")] = py::cast(MERLIN_ALGO_AOBF);
    algorithm_members[py::str("RBFAOO")] = py::cast(MERLIN_ALGO_RBFAOO);
    algorithm_members[py::str("BTE")] = py::cast(MERLIN_ALGO_BTE);
    algorithm_members[py::str("CTE")] = py::cast(MERLIN_ALGO_CTE);

    py::object AlgorithmClass = PyIntEnum("Algorithm", algorithm_members);
    AlgorithmClass.attr("__doc__") = "Probabilistic inference algorithms supported by Merlin.";
    AlgorithmClass.attr("GIBBS").attr("__doc__") = "Gibbs Sampling";
    AlgorithmClass.attr("LBP").attr("__doc__") = "Loopy Belief Propagation";
    AlgorithmClass.attr("IJGP").attr("__doc__") = "Iterative Join Graph Propagation";
    AlgorithmClass.attr("JGLP").attr("__doc__") = "Join Graph Linear Programming";
    AlgorithmClass.attr("WMB").attr("__doc__") = "Weighted Mini-Buckets";
    AlgorithmClass.attr("AOBB").attr("__doc__") = "AND/OR Branch and Bound";
    AlgorithmClass.attr("AOBF").attr("__doc__") = "Best-First AND/OR Search";
    AlgorithmClass.attr("RBFAOO").attr("__doc__") = "Recursive Best-First AND/OR Search";
    AlgorithmClass.attr("BTE").attr("__doc__") = "Bucket-Tree Elimination";
    AlgorithmClass.attr("CTE").attr("__doc__") = "Clique-Tree Elimination";
    m.attr("Algorithm") = AlgorithmClass;

    py::dict task_members;
    task_members[py::str("PR")] = py::cast(MERLIN_TASK_PR);
    task_members[py::str("MAR")] = py::cast(MERLIN_TASK_MAR);
    task_members[py::str("MAP")] = py::cast(MERLIN_TASK_MAP);
    task_members[py::str("MMAP")] = py::cast(MERLIN_TASK_MMAP);
    task_members[py::str("EM")] = py::cast(MERLIN_TASK_EM);

    py::object TaskClass = PyIntEnum("Task", task_members);
    TaskClass.attr("__doc__") = "Probabilistic inference tasks supported by Merlin.";
    TaskClass.attr("PR").attr("__doc__") = "Partition function (probability of evidence)";
    TaskClass.attr("MAR").attr("__doc__") = "Posterior marginals (given evidence)";
    TaskClass.attr("MAP").attr("__doc__") = "Maximum aposteriori (given evidence)";
    TaskClass.attr("MMAP").attr("__doc__") = "Marginal MAP (given evidence)";
    TaskClass.attr("EM").attr("__doc__") = "Parameter learning (EM)";
    m.attr("Task") = TaskClass;

    py::dict input_format_members;
    input_format_members[py::str("MARKOV")] = py::cast(MERLIN_INPUT_MARKOV);
    input_format_members[py::str("BAYES")] = py::cast(MERLIN_INPUT_BAYES);

    py::object InputFormatClass = PyIntEnum("InputFormat", input_format_members);
    InputFormatClass.attr("__doc__") = "Input graphical model formats supported by Merlin.";
    InputFormatClass.attr("MARKOV").attr("__doc__") = "UAI Markov Random Field (default)";
    InputFormatClass.attr("BAYES").attr("__doc__") = "UAI Bayes network";
    m.attr("InputFormat") = InputFormatClass;

    py::dict output_format_members;
    output_format_members[py::str("UAI")] = py::cast(MERLIN_OUTPUT_UAI);
    output_format_members[py::str("JSON")] = py::cast(MERLIN_OUTPUT_JSON);

    py::object OutputFormatClass = PyIntEnum("OutputFormat", output_format_members);
    OutputFormatClass.attr("__doc__") = "Output formats for inference results.";
    OutputFormatClass.attr("UAI").attr("__doc__") = "UAI output format (default)";
    OutputFormatClass.attr("JSON").attr("__doc__") = "JSON output format";
    m.attr("OutputFormat") = OutputFormatClass;

}