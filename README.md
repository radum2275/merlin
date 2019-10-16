# Contact

Radu Marinescu (`radu.marinescu@ie.ibm.com`)

# Description

Merlin is a standalone solver written in `C++` that implements state-of-the-art exact 
and approximate algorithms for probabilistic inference over graphical models 
including both directed and undirected models (e.g., Bayesian networks,
Markov Random Fields). Merlin supports the most common probabilistic inference tasks such as
computing the partition function or probability of evidence (PR), posterior 
marginals (MAR), as well as MAP (also known as maximum aposteriori or most 
probable explanation) and Marginal MAP configurations. In addition, Merlin supports
maximum likelihoood (EM) parameter learning for Bayesian networks only.

Merlin implements the classic Loopy Belief Propagation (LBP) algorithm as well
as more advanced generalized belief propagation algorithms such as Iterative
Join-Graph Propagation (IJGP) and Weighted Mini-Bucket Elimination (WMB). The 
WMB(i) algorithm is parameterized by an i-bound that allows for a controllable 
tradeoff between accuracy of the results and the computational
cost. Larger values of the i-bound typically yield more accurate results but
it takes more time and memory to compute them. Selecting a large enough i-bound
allows for exact inference (i.e., i-bound equal to the treewidth of the model).
For relatively small i-bounds Merlin performs approximate inference.

# Dependencies

Merlin requires a recent version of the `boost` library, and it must be linked
with the `boost_program_options` library.

# API


## Merlin API
Class `Merlin` defined in `merlin.h` header exposes most of the functionality 
of the solver. A graphical model is a collection of factors (or positive
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
* `MERLIN_TASK_EM`    : EM parameter learning for Bayes nets 

        bool set_algorithm(size_t a)
This method sets the the algorithm to be used when solving the selected 
probabilistic inference task. The possible values for the `a` parameter are:
* `MERLIN_ALGO_GIBBS`     : Gibbs sampling
* `MERLIN_ALGO_LBP`       : Loopy belief propagation
* `MERLIN_ALGO_IJGP`      : Iterative join graph propagation
* `MERLIN_ALGO_JGLP`      : Join graph linear programming
* `MERLIN_ALGO_WMB`       : Weighted mini-bucket elimination
* `MERLIN_ALGO_BTE`       : Bucket tree elimination
* `MERLIN_ALGO_CTE`		: Clique Tree Elimination algorithm
* `MERLIN_ALGO_AOBB`      : AND/OR branch and bound search (not implemented yet)
* `MERLIN_ALGO_AOBF`      : Best-first AND/OR search (not implemented yet)
* `MERLIN_ALGO_RBFAOO`    : Recursive best-first AND/OR search (not implemented yet)

        void set_ibound(size_t ibound)
This method sets the i-bound parameter which is used by the following
algorithms: `WMB`, `IJGP`, `JGLP` (as well as search based ones `AOBB`, `AOBF`,
and `RBFAOO`). The default value is `4`.

        void set_iterations(size_t iter)
This method sets the number of iterations to be executed by the inference or learning
algorithm. The parameter is specific to the following algorithms: `WMB`, `IJGP`,
`JGLP`, `LBP` and `EM`. The default value is `100`. For Gibbs sampling consider
runnig several thousands of iterations.

        void set_samples(size_t s)
This method sets the number of samples to be generated in each iteration of the
`GIBBS` sampling algorithm. The default value is `100`. 
 
        void run()
This method runs the inference algorithm for the selected task on the input
graphical model and evidence (if any). The output is generated into a file
specified using the UAI format. The name of the output file is obtained from 
the input file augmented with the `task.out` suffix, where `task` corresponds
to one of the follwing: `PR`, `MAR`, `MAP`, `MMAP` or `EM`. 

 
# Source Code

The source code is organized along the following directory structure and 
requires a standard GNU build using the GNU Autotools toolchain.

* `src` - contains the source files
* `include` - contains the header files
* `data/` - contains several example graphical models
* `doc/` - contains the documentation
* `build/` - contains the intermediate build files
* `bin/` - contains the compiled binary

# Build

The simplest way to compile the solver is to run `make` in the root folder. The 
`Makefile` distributed with this version compiles the solver in `debug` mode
(i.e., `-g` option) on a MacOS platform. For `Ubuntu`, please discard the `-stdlib=libc++`
option in line 2 of the `Makefile`.

## Building the Documentation
Merlin uses Doxygen to build automatically the reference manual of the library,
and supports both `html` and `latex` (see the corresponding `doc/html` and
`doc/latex` subfolders). 

To build the entire documentation, simply run `doxygen merlin.doxygen` in the
main folder `merlin/`. To generate the pdf run `make all` in the
`doc/latex` subfolder).

# Runnig the Solver

Merlin has a very simple command line format and accepts the following command
line arguments:


* `--input-file <filename>` - is the input graphical model file
* `--evidence-file <filename>` - is the evidence file
* `--virtual-evidence-file <filename>` - is the virtual (or likelihood) evidence file
* `--query-file <filename>` - is the query file and is relevant to the MAR and MMAP tasks only
* `--dataset-file <filename>` - is the training dataset file and is relevant to the EM task only
* `--output-file <filename>` - is the output file where the solutions is written
* `--output-format <format>` - is the format of the output file, use either `uai` or `json`
* `--task <task>` - is the inference task and it can be one of the following: `PR`, `MAR`, `MAP`, `MMAP`, `EM`
* `--algorithm <alg>` - is the inference algorithm and it can be one of the following: `wmb`, `bte`, `cte`, `ijgp`,`jglp`, `lbp` or `gibbs`
* `--ibound <n>` - is the ibound used by `wmb`, `ijgp`, `jglp` algorithms
* `--iterations <n>` - is the number of iterations used by `wmb`, `ijgp`, `lbp` and `jglp`
* `--samples <n>` - is the number of samples used by `gibbs`
* `--help` - lists the options
* `--debug` - activate the debug mode
* `--verbose <n>` - is the verbosity level (`1=low`, `2=medium`, `3=high`)
* `--positive` - activate positive mode (i.e., probabilities > 0)
* `--threshold` - is the threshold value used for checking the convergence of an algorithm (default is `1e-06`)
* `--init-factors` - is the factor initialization method (possible values are `none`, `uniform` - default, `random`)
* `--alpha` - is the equivalent sample size used for initializing Dirichlet priors (default is `5.0`)

Run `./merlin --help` to get the complete list of arguments.

Example of command line:

`./merlin --input-file pedigree1.uai --evidence-file pedigree1.evid --task MAR --algorithm wmb --ibound 4 --iterations 10 --output-format json`

This example will run the WMB algorithm with ibound 10 and 10 iterations to solve the MAR task (calculate posterior marginals) on the `pedigree1.uai` instance given evidence `pedigree1.evid`.

All files (input, evidence, query) must be specified in the UAI files format.
See the next section for more details.

# Parameter Learning via EM (Maximum Likelihood Estimation)

For parameter learning we can use the following command line example:

`--input-file cancer.uai --dataset-file cancer.dat --task EM --algorithm cte --iterations 10 --threshold 0.000001`

In this case, Merlin will run EM parameter learning. Note that `cancer.uai` must be a Bayesian network. The file `cancer.dat` contains the training dataset with potentially missing values (see the next section for details on the 
dataset file format).

Merlin performs exact inference by running the Clique Tree Elimination (`CTE`) algorithm.

Merlin allows the user to control the initialization strategy of the factors. Use
`--init-factors` to initialize the factors uniformly (`uniform`) or randomly (`random`).
If the chosen value is `none`, then EM will use the initial factor values provided
in the input model (and assumed to be expert knowledge).

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

## Virtual Evidence File Format
Virtual or likelihood evidence can be specified in a separate file. The first line
contains the number of virtual evidence variable. The subsequent lines correspond
to the virtual evidence variable. For each evidence variable, the line contains
the index of the variable, the domain size and the likelihood values corresponding
to each of the domain values (all numbers specified on a line must be separated
by a single space).

Going back to our example, virtual evidence on variables `Y` and `Z` can be 
specified as follows:

		2
		1 2 0.6 0.8
		2 2 0.1 0.3

## Training Dataset File Format for EM Parameter Learning
The training examples are specified in a CSV file. Each training example
consists of a single comma separated line containing the the domain value index 
of the corresponding variable. Missing values are denoted by `?`.

    	1,0,?,0,?
		1,?,0,?,0

In this example, we have two training examples over 5 variables. Looking at the
first line we see that the value of the first variable is `1`, the value of the 
second variable is `0`, the value of the third variable is missing (`?`), the 
value of the fourth variable is `0`, and the value of the fifth variable is missing (`?`).  

For EM parameter learning, Merlin supports virtual evidence specified in the
training dataset file. Virtual evidence for a variable has to be enclosed between
square brackets `[`, `]` and contains a `;` separated string that specifies the 
likelihood of each domain value.

		1,[0.4;0.8],?,[0.9;0.1],?
    
## Query File Format

The Query file can be used for MAR and Marginal MAP.

### MAR
In this case, the Query file is used to specify a joint marginal. The file consists
of a single line that begins with the number of variables and is followed by the
indices of the query variables.

For example, a query file like:

        2 0 2
specifies a joint marginal over variables `X` and `Z` in our example Markov network.


### Marginal MAP
In this case, the query variables for Marginal MAP inference (i.e., MAP variables)
are specified in a separate file. 
The query file consists of a single line. The line will begin with the number of 
query variables, followed by the indexes of the query variables. The indexes 
correspond to the ones implied by the original problem file.

If, for our example Markov network, we want to use `Y` as the query variable 
the query file would contain the following line:

        1 1

## UAI Output File Format
The first line must contain only the task solved: `PR|MAP|MAR|MMAP`. The rest 
of the file will contain the solution for the task. The format of 
the `<SOLUTION>` part will be described below.

        PR
        <SOLUTION>

The solution format are as follows depending on the task:

### Partition function `PR`
Line with the value of the log (ie, natural logarithm) of the partition function. 
For example, an approximation `ln Pr(e) = -0.2008` which is known to be an upper 
bound may have a solution line:

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

### EM Parameter Learrning `EM`
The output of the EM algorithm is a new model file containing the parameters
learned from the training dataset. The name of the output file is obtained by
appending the `.EM` suffix to the input model filename (e.g., for input filename `cancer.uai`
Merlin generates the output in the file `cancer.uai.EM`).
        
## JSON Output File Format
Merlin supports a JSON format for the output file.

### Partition function `PR`

    { 
        "algorithm" : "wmb",
        "ibound" : 2,  
        "iterations" : 10,  
        "task" : "PR",  
        "value" : -2.551383,  
        "status" : "true",  
    }

### Marginals `MAR`

    { 
        "algorithm" : "lbp",  
        "iterations" : 860,  
        "task" : "MAR",  
        "value" : -2.551383,  
        "status" : "true",  
        "solution" : [ 
            { 
                "variable" : 0,  
                "states" : 2,  
                "marginal" : [0.960694, 0.039306] 
            }, 
            { 
                "variable" : 1,  
                "states" : 2,  
                "marginal" : [0.912524, 0.087476] 
            }
        ]
    }

### Maximum aposteriori `MAP`
The solution contains all variables in the input graphical model (including 
the evidence variables)

    { 
        "algorithm" : "jglp",  
        "ibound" : 2,  
        "iterations" : 10,  
        "task" : "MAP",  
        "value" : -8.801573,  
        "status" : "true",  
        "solution" : [ 
            { 
                "variable" : 0, 
                "value" : 0
            }, 
            { 
                "variable" : 1, 
                "value" : 0
            }
        ] 
    }

### Marginal MAP `MMAP`
The solution contains only the query variables, indexed as in the input query file.

    { 
        "algorithm" : "wmb",  
        "ibound" : 2,  
        "iterations" : 10,  
        "task" : "MMAP",  
        "value" : -12.801573,  
        "status" : "true",  
        "solution" : [ 
            { 
                "variable" : 0, 
                "value" : 0
            }, 
            { 
                "variable" : 1, 
                "value" : 0
            }
        ] 
    }
    
   