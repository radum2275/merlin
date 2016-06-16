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
/// \author Radu Marinescu
///

// Software version
#define VERSIONINFO "libmerlin 1.0.0"
#define COPYRIGHT "(c) Copyright IBM Corp. 2015, 2016\nAll Rights Reserved"

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
#include <malloc.h>
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
#define MERLIN_DOUBLE_PRECISION 6		///<  Precision used for displaying doubles (default 6)


#endif /* IBM_MERLIN_BASE_H_ */
