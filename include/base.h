/*
 * base.h
 *
 *  Created on: Feb 8, 2013
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

/// \file base.h
/// \brief Global definitions
/// \author Radu Marinescu radu.marinescu@ie.ibm.com
///

// Software version
#define VERSIONINFO "libmerlin 1.7.0"
#define COPYRIGHT "(c) Copyright IBM Corp. 2015 - 2019\nAll Rights Reserved"

#ifndef IBM_MERLIN_BASE_H_
#define IBM_MERLIN_BASE_H_

// Debugging purposes
//#define MERLIN_DEBUG

// Standard includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <memory.h>
#include <sys/types.h>
#include <sys/timeb.h>

// STL kernel
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <map>
#include <functional>
#include <algorithm>
#include <deque>
#include <list>
#include <queue>
#include <set>
#include <stack>
#include <exception>
#include <stdexcept>
#include <string>
#include <iomanip>
#include <limits>
#include <numeric>
#include <cmath>

/// Miscelaneous constants
#define MERLIN_PRECISION 	6			///< Precision used for displaying doubles (default 6)
#define MERLIN_EPSILON 		1e-6		///< Small epsilon value to control determinism
#define MERLIN_UNKNOWN		-1			///< Unknown value
#define MERLIN_INIT_RANDOM	10			///< Initialize factors randomly
#define MERLIN_INIT_UNIFORM	11			///< Initialize factors uniformly
#define MERLIN_INIT_NONE	12			///< No factor initialization
#define MERLIN_MAXSIZE_JOINT_MARGINAL 5	///< Maximum size of a joint marginal


///
/// Probabilistic inference algorithms.
///
#define MERLIN_ALGO_GIBBS 	1000		///< Gibbs Sampling
#define MERLIN_ALGO_LBP		1001		///< Loopy Belief Propagation
#define MERLIN_ALGO_IJGP	1002		///< Iterative Join Graph Propagation
#define MERLIN_ALGO_JGLP	1003		///< Join Graph Linear Programming
#define MERLIN_ALGO_WMB		1004		///< Weighted Mini-Buckets
#define MERLIN_ALGO_AOBB	1005		///< AND/OR Branch and Bound
#define MERLIN_ALGO_AOBF	1006		///< Best-First AND/OR Search
#define MERLIN_ALGO_RBFAOO	1007		///< Recursive Best-First AND/OR Search
#define MERLIN_ALGO_BTE		1008		///< Bucket-Tree Elimination
#define MERLIN_ALGO_CTE		1009		///< Clique-Tree Elimination

///
/// Probabilistic inference tasks.
///
#define MERLIN_TASK_PR		10			///< Partition function (probability of evidence)
#define MERLIN_TASK_MAR		20			///< Posterior marginals (given evidence)
#define MERLIN_TASK_MAP		30			///< Maximum aposteriori (given evidence)
#define MERLIN_TASK_MMAP	40			///< Marginal MAP (given evidence)
#define MERLIN_TASK_EM		50			///< Parameter learning (EM)

///
/// Input graphical models.
///
#define MERLIN_INPUT_MARKOV	1			///< UAI Markov Random Filed (default)
#define MERLIN_INPUT_BAYES  2			///< UAI Bayes network

///
/// Output format
///
#define MERLIN_OUTPUT_UAI	10			///< UAI output format (default)
#define MERLIN_OUTPUT_JSON	11			///< JSON output format

#endif /* IBM_MERLIN_BASE_H_ */
