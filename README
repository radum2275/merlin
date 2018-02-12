# Description

Merlin is an extensible `C++` library that implements state-of-the-art exact 
and approximate algorithms for probabilistic inference over graphical models 
including both directed and undirected models (e.g., Bayesian networks,
Markov Random Fields). It can be used for many applications and research in
bioinformatics, computer vision, or speech and language processing to name a
few. Merlin supports the most common probabilistic inference tasks such as
computing the partition function or probability of evidence (PR), posterior 
marginals (MAP), as well as MAP (also known as maximum aposteriori or most 
probable explanation) and marginal MAP configurations.

Merlin implements the classic Loopy Belief Propagation (LBP) algorithm as well
as more advanced generalized belief propagation algorithms such as Iterative
Join-Graph Propagation (IJGP) and Weighted Mini-Bucket Elimination (WMB). The
default algorithm for computing all four probabilistic inference tasks (PR, MAR,
MAP, MMAP) is WMB(i). The algorithm is parameterized by an i-bound that allows
for a controllable tradeoff between accuracy of the results and the computational
cost. Larger values of the i-bound typically yield more accurate results but
it takes more time and memory to compute them. Selecting a large enough i-bound
allows for exact inference (i.e., i-bound equal to the treewidth of the model).
For relatively small i-bounds Merlin performs approximate inference.


# API

## Basic (Black-Box) API
Merlin exposes (see `merlin.h` header) a single C-like function called `run()`:

        int run(const char* inputFile, 
            const char* evidenceFile, 
            const char* queryFile, 
            const char* outputFile, 
            const char* task, 
            const unsigned int ibound, 
            const unsigned int iterations);

where:
* `inputFile` - name of the file containing the input graphical model specified in the UAI format (see main documentation for details on the format).
* `evidenceFile` - name of the file containing the evidence or observations as variable and value pairs (variable indeces are the same as in the `inputFile`).
* `queryFile` - name of the file containing the query variable for the marginal MAP inference task (variable indeces are the same as in the `inputFile`).
* `outputFile` - name of the file containing the result of the inference task (see main documentation for details on the format).
* `task` - name of the probabilistic inference task. The following names can be used: PR, MAR, MAP, MMAP.
* `ibound` - parameter i-bound of the inference algorithm `WMB(i)` controling the accuracy of the results. The value should be greater or equal to 2 (default is 2).
* `iterations` - parameter iterations of the inference algorithm `WMB(i)` specifing the number of iterations to be executed (default is 100).

The function returns 0 if successful, in which case the output file contains the
solution to the probabilistic inference task. Otherwise, the return value is 1, 
in which case the content of the output file is undefined and it should be
ignored.


### Example

A simple application that uses the library is included with the distribution 
(see the `demo.cpp` file in the `example/` folder). The source code is given 
below. The input graphical model is given in the `pedigree1.uai` file, while 
`pedigree1.evid` contains a set of observations (ie, variable value pairs). 
In addition, the file `pedigree1.map` contains the query variables for the 
marginal MAP task. The application solves the following inference tasks (the 
results are written in the corresponding `.out` files):

* MAR = posterior marginals (ie, marginal distributions of all unobserved variables given the evidence)
* MAP = maximum aposteriori (ie, configuration of unobserved variables having maximum probability)
* MMAP = marginal MAP (ie, configuration of the query variables having maximum marginal probability)

        #include "merlin.h"
    
        void demo_run(void) {
        
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

The application should be compiled with the `-lmerlin` option. It is also
possible to just run the `example` wrapper script (generated automatically
during by the build scripts) in the `example/` folder.


## Advanced API
Class `Merlin` defined in `merlin.h` header exposes most of the functionality 
of the library. A graphical model is a collection of factors (or positive
real-valued functions) defined over subsets of variables. Variables are
assumed to be indexed from `0`.

### Methods

        bool read_model(const char* f)
This method loads the graphical model from a file which is specified using the
UAI format (see also the File Formats section). Returns `true` if successful
and `false` otherwise.

        bool read_evidence(const char* f)
This method loads the evidence variables and their corresponding observed values
from a file which is also specified using the UAI format. Returns `true` if
successful, and `false` otherwise.

        bool read_query(const char* f)
This method loads the query variables from a file specified using the UAI format.
The query variables (also known as MAX of MAP variables) are only specific to 
Marginal MAP (MMAP) inference tasks. Returns `true` if successful, and `false` 
otherwise.

        bool set_task(size_t t)
This method sets the probabilistic inference task to be solved. The possible
values for the `t` parameter are:
* `MERLIN_TASK_PR`    : Partition function (probability of evidence)
* `MERLIN_TASK_MAR`   : Posterior marginals (given evidence)
* `MERLIN_TASK_MAP`   : Maximum aposteriori (given evidence)
* `MERLIN_TASK_MMAP`  : Marginal MAP (given evidence)

        bool set_algorithm(size_t a)
This method sets the the algorithm to be used when solving the selected 
probabilistic inference task. The possible values for the `a` parameter are:
* `MERLIN_ALGO_GIBBS`     : Gibbs sampling
* `MERLIN_ALGO_LBP`       : Loopy belief propagation
* `MERLIN_ALGO_IJGP`      : Iterative join graph propagation
* `MERLIN_ALGO_JGLP`      : Join graph linear programming
* `MERLIN_ALGO_WMB`       : Weighted mini-bucket elimination
* `MERLIN_ALGO_AOBB`      : AND/OR branch and bound search (not implemented)
* `MERLIN_ALGO_AOBF`      : Best-first AND/OR search (not implemented)
* `MERLIN_ALGO_RBFAOO`    : Recursive best-first AND/OR search (not implemented)

        void set_param_ibound(size_t ibound)
This method sets the i-bound parameter which is used by the following
algorithms: `WMB`, `IJGP`, `JGLP` (as well as search based ones `AOBB`, `AOBF`,
and `RBFAOO`). The default value is `4`.

        void set_param_iterations(size_t iter)
This method sets the number of iterations to be executed by the inference 
algorithm. The parameter is specific to the following algorithms: `WMB`, `IJGP`,
`JGLP`, and `LBP`. The default value is `100`. For Gibbs sampling consider
runnig several thousands of iterations.

        void set_param_samples(size_t s)
This method sets the number of samples to be generated in each iteration of the
`GIBBS` sampling algorithm. The default value is `100`. 
 
        void run()
This method runs the inference algorithm for the selected task on the input
graphical model and evidence (if any). The output is generated into a file
specified using the UAI format. The name of the output file is obtained from 
the input file augmented with the `task.out` suffix, where `task` corresponds
to one of the follwing: `PR`, `MAR`, `MAP`, or `MMAP`. 

### Example
A simple example using the advanced Merlin API is included in the `example/`folder
(i.e., `demo.cpp` file). The source code is given below. The input graphical model 
is given in the `pedigree1.uai` file, while `pedigree1.evid` contains a set 
of observations (ie, variable value pairs). In addition, the file `pedigree1.map` 
contains the query variables for the Marginal MAP task. The application solves all four
inference tasks and the results are written in the corresponding `.out` files.
    
        #include "merlin.h"
        
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

The application is automatically compiled with the `-lmerlin` option. It is also
possible to just run the `example` wrapper script (generated automatically
during by the build scripts) in the `example/` folder.

 
# Source Code

The source code is organized along the following directory structure and 
requires a standard GNU build using the GNU Autotools toolchain.

* `src/` - contains the source (.cpp) files
* `include/` - contains the header (.h) files
* `example/` - contains an example program (demo.cpp) that uses the library
* `doc/` - contains the documentation

# Build

The simplest way to compile the library is to `cd` to the directory containing 
the source code and type `./configure`. Then type `make` to compile the library.
The binaries will be generated in the default location, `src/.libs/`. `make` will
also build the example program in `demo/`. Finally, type `make install` to install the
library in the default location (`/usr/lib` for the `libmerlin.so` object and `/usr/include` 
for `merlin.h`). When installing into a prefix owned by the root, it is recommended 
that the library be configured and built as a regular user, and only the 
`make install` phase executed with root privileges (see INSTALL for more details).

        -$ ./configure
        -$ make
        -$ sudo make install
    
The interface of the library is exposed in the `merlin.h` header file which must be
included in the source files of the application. See demo/ for a simple example.
The name of the shared library is `libmerlin`.

## Building the Documentation
Merlin uses Doxygen to build automatically the reference manual of the library,
and supports both `html` and `latex` (see the corresponding `doc/html` and
`doc/latex` subfolders). 

To build the entire documentation, simply run `doxygen merlin.doxygen` in the
main folder `merlin/`. To generate the pdf run `make all` in the
`doc/latex` subfolder).


# File Formats

## Input File Format

Merlin uses a simple text file format which is specified below to describe a 
problem instances (i.e., graphical model). The format is identical to the one 
used during the UAI Inference competitions.

### Structure
The input file format consists of the following two parts, in that order:

        <Preamble>
        <Factors>

The contents of each section (denoted <…> above) are described in the following:

#### Preamble
The description of the format will follow a simple Markov network with three 
variables and two functions. A sample preamble for such a network is:

        MARKOV
        3
        2 2 3
        2
        2 0 1
        2 1 2

The preamble starts with one line denoting the type of network. Generally, this 
can be either BAYES (if the network is a Bayesian network) or MARKOV (in case of 
a Markov network).

The second line contains the number of variables.

The third line specifies the cardinalities of each variable, one at a time, 
separated by a whitespace (note that this implies an order on the variables 
which will be used throughout the file).

The fourth line contains only one integer, denoting the number of cliques in the 
problem. Then, one clique per line, the scope of each clique is given as follows: 
The first integer in each line specifies the number of variables in the clique, 
followed by the actual indexes of the variables. The order of this list is not 
restricted (for Bayesian networks we assume that the child variable of the clique 
is the last one). Note that the ordering of variables within a factor will 
follow the order provided here.

Referring to the example above, the first line denotes the Markov network, the 
second line tells us the problem consists of three variables, let's refer to 
them as `X`, `Y`, and `Z`. Their cardinalities are `2`, `2`, and `3` 
respectively (from the third line). Line four specifies that there are 2 cliques. 
The first clique is `X,Y`, while the second clique is `Y,Z`. Note that 
variables are indexed starting with `0`.

#### Factors
Each factor is specified by giving its full table (i.e, specifying a 
non-negative real value for each assignment). The order of the factors is 
identical to the one in which they were introduced in the preamble, the first 
variable has the role of the 'most significant' digit. For each factor table, 
first the number of entries is given (this should be equal to the product of the 
domain sizes of the variables in the scope). Then, one by one, separated by 
whitespace, the values for each assignment to the variables in the factor's 
scope are enumerated. Tuples are implicitly assumed in ascending order, with 
the last variable in the scope as the `least significant`. To illustrate, we 
continue with our Markov network example from above, let's assume the following 
conditional probability tables:

        X | P(X)  
        0 | 0.436 
        1 | 0.564 
        
        X   Y |  P(Y,X)
        0   0 |  0.128
        0   1 |  0.872
        1   0 |  0.920
        1   1 |  0.080
        
        Y   Z |  P(Z,Y)
        0   0 |  0.210
        0   1 |  0.333
        0   2 |  0.457
        1   0 |  0.811
        1   1 |  0.000
        1   2 |  0.189

Then we have the corresponding file content:

        2
         0.436 0.564
        
        4
         0.128 0.872
         0.920 0.080
        
        6
         0.210 0.333 0.457
         0.811 0.000 0.189

Note that line breaks and empty lines are effectively just a whitespace, 
exactly like plain spaces “ ”. They are used here to improve readability.


## Evidence File Format
Evidence is specified in a separate file. The evidence file consists of a single
line. The line will begin with the number of observed variables in the sample, 
followed by pairs of variable and its observed value. The indexes correspond to 
the ones implied by the original problem file.

If, for our example Markov network, `Y` has been observed as having its first 
value and `Z` with its second value, the evidence file would contain the 
following line:

        2 1 0 2 1
    
## Query File Format
Query variables for Marginal MAP inference are specified in a separate file. 
The query file consists of a single line. The line will begin with the number of 
query variables, followed by the indexes of the query variables. The indexes 
correspond to the ones implied by the original problem file.

If, for our example Markov network, we want to use `Y` as the query variable 
the query file would contain the following line:

        1 1

## Output File Format
The first line must contain only the task solved: `PR|MAP|MAR|MMAP`. The rest 
of the file will contain the solution for the task. Solvers can write more then 
one solution by writing `-BEGIN-` at the head of the new solution. In the 
example below the task we choose is `PR`. We have two solutions. The format of 
the `<SOLUTION>` part will be described below.

        PR
        <SOLUTION>
        -BEGIN-
        <SOLUTION>

The solution format are as follows depending on the task:

### Partition function `PR`
Line with the value of the log10 of the partition function. For example, an 
approximation `log10 Pr(e) = -0.2008`which is known to be an upper bound may 
have a solution line:

        -0.2008
### Maximum aposteriori `MAP`
A space separated line that includes:
* the number `n` of model variables, and
* the MAP instantiation, a list of value indices for all `n` variables.

For example, an input model with 3 binary variables may have a solution line:

        3 0 1 0

### Marginals `MAR`
A space separated line that includes:
* the number of variables `n` in the model, and
* a list of marginal approximations of all the variables. For each variable 
its cardinality is first stated, then the probability of each state is stated. 
The order of the variables is the same as in the model, all data is space 
separated.

For example, a model with `3` variables, with cardinalities of `2`, `2`, `3`
respectively. The solution might look like this:

        3 2 0.1 0.9 2 0.3 0.7 3 0.2 0.2 0.6

### Marginal MAP `MMAP`
A space separated line that includes:
* the number `q` of query (or MAP) variables, and
* their most probable instantiation, a list of variable value pairs for all `q`
variables.

For example, if the solution is an assignment of `0`, `1` and `0` to three query 
variables indexed by `2` `3` and `4` respectively, the solution will look as follows:

        3 2 0 3 1 4 0
